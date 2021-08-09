












! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ocn_rpn_calculator
!
!> \brief MPAS ocean analysis core member: rpn_calculator
!> \author Jon Woodring
!> \date   March 21, 2016
!> \details
!>  Flexible vector RPN calculator of MPAS fields for up to 2D fields.
!-----------------------------------------------------------------------
module ocn_rpn_calculator
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  use ocn_constants

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: &
         ocn_init_rpn_calculator, &
         ocn_compute_rpn_calculator, &
         ocn_restart_rpn_calculator, &
         ocn_finalize_rpn_calculator

  ! Private module variables
  !--------------------------------------------------------------------

  type rpn_stack_value_type
    integer :: symbol_type
    integer :: number_of_dims

    type (field0DReal), pointer :: d0
    type (field1DReal), pointer :: d1
    type (field2DReal), pointer :: d2
  end type rpn_stack_value_type

  integer, parameter :: SYMBOL_NOT_FOUND = 0

  integer, parameter :: IS_OPERATOR = 10
  integer, parameter :: IS_VARIABLE = 100
  integer, parameter :: IS_TEMPORARY = 1000

  integer, parameter :: MAX_STACK_SIZE = StrKIND / 2

  character (len=1), dimension(8), parameter :: variable_names = &
    (/ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' /)
  character (len=3), dimension(4) :: operator_names = &
    (/ '*  ' , '+  ', '-  ', '/  ' /)
  integer, parameter :: MUL_OP = 1
  integer, parameter :: PLUS_OP = 2
  integer, parameter :: MINUS_OP = 3
  integer, parameter :: DIV_OP = 4
  ! TODO FIXME
  ! integer, parameter :: SUM_OP = 5

  character (len=1), dimension(4) :: expression_names = &
    (/ '1', '2', '3', '4' /)

  character (len=StrKIND), parameter :: VARIABLE_PREFIX = &
    'config_AM_rpnCalculator_variable_'
  character (len=StrKIND), parameter :: EXPRESSION_PREFIX = &
    'config_AM_rpnCalculator_expression_'
  character (len=StrKIND), parameter :: OUTPUT_PREFIX = &
    'config_AM_rpnCalculator_output_name_'

  character (len=StrKIND), parameter :: OUTPUT_STREAM_CONFIG = &
    'config_AM_rpnCalculator_output_stream'

  character (len=StrKIND), parameter :: NONE_TOKEN = 'none'

  character (len=StrKIND), parameter :: MPAS_CORE_NAME = 'MPAS-Ocean'

!***********************************************************************
contains



!***********************************************************************
! routine ocn_init_rpn_calculator
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, last, stack_pointer
  character (len=StrKIND) :: config, field_name
  character (len=StrKIND), pointer :: config_result
  type (rpn_stack_value_type) :: output_value
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE) :: stack

  ! start procedure
  err = 0

  ! typecheck all the expressions
  last = size(expression_names)
  do i = 1, last
    config = trim(EXPRESSION_PREFIX) // trim(expression_names(i))
    call mpas_pool_get_config(domain % configs, config, config_result)

    if (trim(config_result) /= trim(NONE_TOKEN)) then
      stack_pointer = -1 ! typecheck with an empty stack
      call eval_expression(domain, config_result, i, stack, stack_pointer)

      ! check the stack size
      if (stack_pointer /= 1) then
        call mpas_log_write( &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM ' // &
          'resulted in the stack size not being equal to 1: ' // &
          'i.e., the return result of the expression should be the only ' // &
          'value on the stack after evaluation', MPAS_LOG_CRIT)
      end if

      ! check that it's a new value
      if (stack(stack_pointer) % symbol_type /= IS_TEMPORARY) then
        call mpas_log_write( &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM did not calculate anything, ' // &
          ' i.e., it only pushed a variable onto the stack', MPAS_LOG_CRIT)
      end if

      ! rename the stack field and put in allFields pool
      config = trim(OUTPUT_PREFIX) // trim(expression_names(i))
      call mpas_pool_get_config(domain % configs, config, config_result)

      if (trim(config_result) == (NONE_TOKEN)) then
        call mpas_log_write( &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM was set, but the output field name ' // &
          'for that expression was set to "none"', MPAS_LOG_CRIT)
      end if

      if (stack(1) % number_of_dims == 0) then
        stack(1) % d0 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d0)
      else if (stack(1) % number_of_dims == 1) then
        stack(1) % d1 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d1)
      else if (stack(1) % number_of_dims == 2) then
        stack(1) % d2 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d2)
      else
        call mpas_log_write( &
          'the impossible happened, the dimensions of the result on the ' // &
          'stack, for expression #' // trim(expression_names(i)) // &
          ' was not between 0 and 2 in the RPN calculator AM', MPAS_LOG_CRIT)
      end if
      field_name = config_result

      ! put them in the stream if necessary
      call mpas_pool_get_config(domain % configs, &
        OUTPUT_STREAM_CONFIG, config_result)

      if (trim(config_result) /= trim(NONE_TOKEN)) then
        call mpas_stream_mgr_add_field(domain % streamManager, &
          config_result, field_name, ierr=err)
      end if
    end if
  end do

end subroutine ocn_init_rpn_calculator!}}}



!***********************************************************************
! routine ocn_compute_rpn_calculator
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_rpn_calculator(domain, timeLevel, err)!{{{
  ! input variables
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, stack_pointer, last
  character (len=StrKIND) :: config
  character (len=StrKIND), pointer :: config_result
  type (rpn_stack_value_type) :: output_value
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE) :: stack
  type (field0DReal), pointer :: d0
  type (field1DReal), pointer :: d1, t1
  type (field2DReal), pointer :: d2, t2
  real (kind=RKIND), dimension(:), pointer :: s1
  real (kind=RKIND), dimension(:,:), pointer :: s2

  ! start procedure
  err = 0

  ! do all the expressions
  last = size(expression_names)
  do i = 1, last
    config = trim(EXPRESSION_PREFIX) // trim(expression_names(i))
    call mpas_pool_get_config(domain % configs, config, config_result)

    if (trim(config_result) /= trim(NONE_TOKEN)) then
      stack_pointer = 0 ! evaluate with an empty stack
      call eval_expression(domain, config_result, i, stack, stack_pointer)

      ! lookup the field and reassign pointers - then deallocate stack
      config = trim(OUTPUT_PREFIX) // trim(expression_names(i))
      call mpas_pool_get_config(domain % configs, config, config_result)
      if (stack(1) % number_of_dims == 0) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d0, 1)
        d0 % scalar = stack(1) % d0 % scalar
        call mpas_deallocate_field(stack(1) % d0)
      else if (stack(1) % number_of_dims == 1) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d1, 1)
        t1 => stack(1) % d1
        do while (associated(d1))
          s1 => d1 % array
          d1 % array => t1 % array
          t1 % array => s1

          d1 => d1 % next
          t1 => t1 % next
        end do
        call mpas_deallocate_field(stack(1) % d1)
      else
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d2, 1)
        t2 => stack(1) % d2
        do while (associated(d2))
          s2 => d2 % array
          d2 % array => t2 % array
          t2 % array => s2

          d2 => d2 % next
          t2 => t2 % next
        end do
        call mpas_deallocate_field(stack(1) % d2)
      end if

    end if
  end do

end subroutine ocn_compute_rpn_calculator!}}}



!***********************************************************************
! routine ocn_restart_rpn_calculator
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_restart_rpn_calculator!}}}



!***********************************************************************
! routine ocn_finalize_rpn_calculator
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_finalize_rpn_calculator!}}}

!
! local subroutines
!

!***********************************************************************
! routine eval_expression
!
!> \brief Given a character string, evaluate the stack expression
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a character string, evaluate the stack expression
!> and copy the bottom (top) of a 1-length stack into the target MPAS field.
!-----------------------------------------------------------------------
subroutine eval_expression (domain, expression, exp_number, &
    stack, stack_pointer)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: expression
  integer, intent(in) :: exp_number

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables
  integer :: symbol_type
  logical :: eol, typechecking
  character (len=StrKIND) :: symbol, remainder

  ! start procedure
  if(stack_pointer < 0) then
    typechecking = .true.
    stack_pointer = -1 - stack_pointer
  else
    typechecking = .false.
  end if

  eol = .false.
  remainder = expression

  ! get the first symbol
  call stack_token(symbol, remainder, eol)

  ! iterate over symbols
  do while(.not. eol)
    symbol_type = symbol_table(symbol)

    ! operator
    if ((symbol_type > IS_OPERATOR) .and. (symbol_type < IS_VARIABLE)) then
      call eval_operator(exp_number, &
        symbol_type - IS_OPERATOR, stack, stack_pointer, typechecking)
    ! variable
    else &
    if ((symbol_type > IS_VARIABLE) .and. (symbol_type < IS_TEMPORARY)) then
      call eval_variable(domain, exp_number, &
        symbol_type - IS_VARIABLE, stack, stack_pointer, typechecking)
    ! symbol not found
    else
      call mpas_log_write( &
        trim(symbol) // '" found in expression #' // &
        trim(expression_names(exp_number)) // &
        ' in the RPN calculator AM was not found.', MPAS_LOG_CRIT)
    end if

    ! get the next symbol
    call stack_token(symbol, remainder, eol)
  end do

end subroutine eval_expression!}}}



!***********************************************************************
! routine stack_token
!
!> \brief Get the next stack token given a character string
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Parses a character string to get the next stack token to eval.
!-----------------------------------------------------------------------
subroutine stack_token(substr, next, eol)!{{{
  ! input variables

  ! input/output variables
  character (len=StrKIND), intent(inout) :: next

  ! output variables
  character (len=StrKIND), intent(out) :: substr
  logical, intent(out) :: eol

  ! local variables
  integer :: i
  character (len=StrKIND) :: copy

  ! make a copy
  copy = trim(next)

  ! if there's anything in it other than whitespace, pass through
  i = verify(copy, ' ')
  eol = i < 1
  if (eol) then
    return
  end if
  copy = trim(next(i:))

  ! find the first whitespace and split
  i = scan(copy, ' ')

  ! return that substring and the remainder
  if (i > 0) then
    substr = trim(copy(1:i-1))
    next = trim(copy(i+1:))
  else
    substr = trim(copy)
    next = ''
  end if

end subroutine stack_token!}}}



!***********************************************************************
! function symbol_table
!
!> \brief Tries to find the symbol in the symbol table and its value
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Will attempt to find a symbol in the symbol table. The value
!> results are dependent on the return code of the symbol table lookup.
!-----------------------------------------------------------------------
integer function symbol_table (symbol)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: symbol

  ! input/output variables

  ! local variables
  integer :: i, last

  ! start procedure

  ! check the operations
  last = size(variable_names)
  do i = 1, last
    if (trim(symbol) == trim(variable_names(i))) then
      symbol_table = IS_VARIABLE + i
      return
    end if
  end do

  ! check the variables
  last = size(operator_names)
  do i = 1, last
    if (trim(symbol) == trim(operator_names(i))) then
      symbol_table = IS_OPERATOR + i
      return
    end if
  end do

  ! else not found
  symbol_table = SYMBOL_NOT_FOUND

end function symbol_table!}}}



!***********************************************************************
! routine eval_operator
!
!> \brief Given a operator index number, put it on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a operator index number, put the result on the top of
!> the stack. It will combine whatever is on the stack to be able to
!> generate results and push them into the stack.
!-----------------------------------------------------------------------
subroutine eval_operator (exp_number, &
    op_index, stack, stack_pointer, type_checking)!{{{
  ! input variables
  integer, intent(in) :: exp_number, op_index
  logical, intent(in) :: type_checking

  ! input/output variables
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables

  ! start procedure
  if (op_index == MUL_OP) then
    call mul_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == PLUS_OP) then
    call plus_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == MINUS_OP) then
    call minus_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == DIV_OP) then
    call div_operator(exp_number, stack, stack_pointer, type_checking)
  ! TODO FIXME
  ! sum (and other reduces) needs to be fixed,
  ! because it is using (:) over decomposed dimensions, which is wrong
  !
  ! else if (op_index == SUM_OP) then
  !  call sum_operator(exp_number, stack, stack_pointer, type_checking)
  else
    call mpas_log_write( &
      'the impossible happened, tried to apply an unknown operator ' // &
      'in expression #' // trim(expression_names(exp_number)) // &
      ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if

end subroutine eval_operator!}}}



!***********************************************************************
! routine eval_variable
!
!> \brief Given a variable index number, put it on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a variable index number, put the result on the top of
!> the stack. This will look up the field names in the variable and look
!> it up from the framework to push the pointer onto the stack.
!-----------------------------------------------------------------------
subroutine eval_variable (domain, exp_number, &
    var_index, stack, stack_pointer, type_checking)!{{{
  ! input variables
  integer, intent(in) :: exp_number, var_index
  logical, intent(in) :: type_checking

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables
  character (len=StrKIND) :: config
  character (len=StrKIND), pointer :: config_result
  type (mpas_pool_field_info_type) :: info

  ! start procedure
  config = trim(VARIABLE_PREFIX) // trim(variable_names(var_index))
  call mpas_pool_get_config(domain % configs, config, config_result)

  if (type_checking) then
    if (trim(config_result) == trim(NONE_TOKEN)) then
      call mpas_log_write( &
        'the MPAS field assigned to variable ' // &
        trim(variable_names(var_index)) // ' was evaluated, but it is ' // &
        'currently set to "none"', MPAS_LOG_CRIT)
    end if
  end if

  call mpas_pool_get_field_info &
    (domain % blocklist % allFields, config_result, info)

  ! check if it's real
  if (type_checking) then
    if (info % fieldType /= MPAS_POOL_REAL) then
      call mpas_log_write( &
        'the MPAS field "' // trim(config_result) // &
        '"assigned to variable ' // &
        trim(variable_names(var_index)) // ' in the RPN calculator AM is ' // &
        'not a real field', MPAS_LOG_CRIT)
    end if

    ! check if it's 0D-2D
    if (info % nDims > 2) then
      call mpas_log_write( &
        'the MPAS field "' // trim(config_result) // &
        '"assigned to variable ' // &
        trim(variable_names(var_index)) // ' in the RPN calculator AM is ' // &
        'not a 0D, 1D, or 2D field', MPAS_LOG_CRIT)
    end if
  end if

  ! increment the stack and put it on the stack
  stack_pointer = stack_pointer + 1
  stack(stack_pointer) % number_of_dims = info % nDims
  stack(stack_pointer) % symbol_type = IS_VARIABLE

  ! get the dimension name if it is 1D
  if (info % nDims == 0) then
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d0, 1)
  else if (info % nDims == 1) then
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d1, 1)
  else
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d2, 1)
  end if
end subroutine eval_variable!}}}



!***********************************************************************
! routine create_2d_field_from_1ds
!
!> \brief Generates a new 2D field from 1D fields
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take two 1D fields (second and top) and
!> generate a new 2D field (head) with second's dimension as its
!> first dimension, and top's as it's second dimension. If top is
!> decomposed, head will be decomposed as well. If second has
!> constituent names, head will have constituent names as well.
!> Both fields need to be active, otherwise head will be inactive.
!-----------------------------------------------------------------------
subroutine create_2d_field_from_1ds(second, top_head, head)!{{{
  type (field1DReal), pointer, intent(in) :: second, top_head
  type (field2DReal), pointer, intent(out) :: head

  type (field2DReal), pointer :: dst, prev
  type (field1DReal), pointer :: top

  if (mpas_threading_get_thread_num() == 0 ) then
  nullify(head)
  nullify(prev)

  top => top_head
  do while (associated(top))

    ! allocate the linked list for the field blocks
    allocate(dst) 
    if (.not. associated(head)) then
      head => dst
    end if

    if (.not. associated(prev)) then
      nullify(dst % prev)
    else
      prev % next => dst
      dst % prev => prev
    end if
    nullify(dst % next)

    ! copy field info
    dst % fieldName = trim(second % fieldName) // '_' // trim(top % fieldName)
    dst % isDecomposed = top % isDecomposed

    dst % block => top % block
    dst % isVarArray = second % isVarArray
    dst % defaultValue = second % defaultValue
    dst % isActive = top % isActive .and. second % isActive
    dst % hasTimeDimension = &
      top % hasTimeDimension .or. second % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList
    dst % isPersistent = top % isPersistent .or. second % isPersistent


    ! copy constitutent names if second has them
    if (associated(second % constituentNames)) then
      allocate(dst % constituentNames(size(second % constituentNames, dim=1)))
      allocate(dst % attLists(size(second % constituentNames, dim=1)))

      dst % constituentNames(:) = second % constituentNames(:)
    else
      nullify(dst % constituentNames)
      allocate(dst % attLists(1))
    end if

    dst % dimNames(1) = second % dimNames(1)
    dst % dimNames(2) = top % dimNames(1)
    dst % dimSizes(1) = second % dimSizes(1)
    dst % dimSizes(2) = top % dimSizes(1)

    ! allocate memory
    if (top % isActive .and. second % isActive) then
      allocate(dst % array(size(second % array), size(top % array)))
    else
      nullify(dst % array)
    end if

    top => top % next
    prev => dst
  end do
  end if
end subroutine create_2d_field_from_1ds!}}}



!***********************************************************************
! routine create_1d_field_from_2d
!
!> \brief Generates a new 1D field from a 2D
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take a 2D field (top_head) and
!> generate a new 1D field (head) with top's dimension as its
!> dimension. If top is decomposed, head will be decomposed as well.
!-----------------------------------------------------------------------
subroutine create_1d_field_from_2d(top_head, head)!{{{
  type (field2DReal), pointer, intent(in) :: top_head
  type (field1DReal), pointer, intent(out) :: head

  type (field1DReal), pointer :: dst, prev
  type (field2DReal), pointer :: top

  if (mpas_threading_get_thread_num() == 0 ) then
  nullify(head)
  nullify(prev)

  top => top_head
  do while (associated(top))

    ! allocate the linked list for the field blocks
    allocate(dst) 
    if (.not. associated(head)) then
      head => dst
    end if

    if (.not. associated(prev)) then
      nullify(dst % prev)
    else
      prev % next => dst
      dst % prev => prev
    end if
    nullify(dst % next)

    ! copy field info
    dst % fieldName = '_' // trim(top % fieldName)
    dst % isDecomposed = top % isDecomposed

    dst % block => top % block
    dst % isVarArray = .false.
    dst % defaultValue = top % defaultValue
    dst % isActive = top % isActive 
    dst % hasTimeDimension = top % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList
    dst % isPersistent = top % isPersistent

    allocate(dst % attLists(1))

    nullify(dst % constituentNames)

    dst % dimNames(1) = top % dimNames(2)
    dst % dimSizes(1) = top % dimSizes(2)

    ! allocate memory
    if (top % isActive) then
      allocate(dst % array(size(top % array, 2)))
    else
      nullify(dst % array)
    end if

    top => top % next
    prev => dst
  end do
  end if
end subroutine create_1d_field_from_2d!}}}



!***********************************************************************
! routine create_0d_field_from_1d
!
!> \brief Generates a new 1D field from a 2D
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take a 1D field (top) and
!> generate a new 0D field (dst).
!-----------------------------------------------------------------------
subroutine create_0d_field_from_1d(top, dst)!{{{
  type (field0DReal), pointer, intent(out) :: dst
  type (field1DReal), pointer, intent(inout) :: top

  if (mpas_threading_get_thread_num() == 0 ) then
    ! allocate the linked list for the field blocks
    allocate(dst) 
    nullify(dst % prev)
    nullify(dst % next)

    ! copy field info
    dst % fieldName = '_' // trim(top % fieldName)
    dst % isDecomposed = .false.

    dst % block => top % block
    dst % isVarArray = .false.
    dst % defaultValue = top % defaultValue
    dst % isActive = top % isActive 
    dst % hasTimeDimension = top % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList

    allocate(dst % attLists(1))

    nullify(dst % constituentNames)
  end if
end subroutine create_0d_field_from_1d!}}}



!***********************************************************************
! routine mul_operator
!
!> \brief Do mul on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> multiply them together, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine mul_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d
  op_name = '*'
! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were less than two operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer - 1) % number_of_dims
  b_nd = stack(stack_pointer) % number_of_dims

  ! call the right one
  if (a_nd == 0) then
    if (b_nd == 0) then

! 0d 1d
  call mul_op_0d_0d(stack, stack_pointer)
! 0d 1d

    else if (b_nd == 1) then

! 0d 2d
  call mul_op_0d_1d(stack, stack_pointer)
! 0d 2d

    else

! 1d 0d
  call mul_op_0d_2d(stack, stack_pointer)
! 1d 0d

    end if
  else if (a_nd == 1) then
    if (b_nd == 0) then

! 1d 1d same
  call mul_op_1d_0d(stack, stack_pointer)
! 1d 1d same

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 1d 1d diff
  call mul_op_1d_1d_same(stack, stack_pointer)
! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first
  call mul_op_1d_1d_diff(stack, stack_pointer)
! 1d 2d first

      end if
    else
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d2 % dimNames(1))) then

! 1d 2d second
  call mul_op_1d_2d_first(stack, stack_pointer)
! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d
  call mul_op_1d_2d_second(stack, stack_pointer)
! 2d 0d

      end if
    end if
  else 
    if (b_nd == 0) then

! 2d 1d first
  call mul_op_2d_0d(stack, stack_pointer)
! 2d 1d first

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d2 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 2d 1d second
  call mul_op_2d_1d_first(stack, stack_pointer)
! 2d 1d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
            trim(stack(stack_pointer) % d1 % dimNames(1))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 2d
  call mul_op_2d_1d_second(stack, stack_pointer)
! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end
  call mul_op_2d_2d(stack, stack_pointer)
! end

    end if
  end if
end subroutine mul_operator!}}}

subroutine mul_op_0d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d0, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d0

  ! do operation
  top => top_iter % scalar
  second => second_iter % scalar
  temp_iter % scalar = &

! 1-2 break
  second * top
! 1-2 break

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
end subroutine mul_op_0d_0d!}}}

subroutine mul_op_0d_1d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d0

  second => second_iter % scalar

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine mul_op_0d_1d!}}}

subroutine mul_op_0d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d0

  second => second_iter % scalar

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_0d_2d!}}}

subroutine mul_op_1d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d1

  top => top_iter % scalar 

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array 
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine mul_op_1d_0d!}}}

subroutine mul_op_1d_1d_same (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine mul_op_1d_1d_same!}}}

subroutine mul_op_1d_1d_diff (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, j, iend, jend

  ! allocate a temp for result
  call create_2d_field_from_1ds( &
    stack(stack_pointer - 1) % d1, stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array
  iend = size(second)

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top)

    do j = 1, jend
      do i = 1, iend
        temp_iter % array(i,j) = &

! 1-2 break
  second(i) * top(j)
! 1-2 break

      end do
    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_1d_1d_diff!}}}

subroutine mul_op_1d_2d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  second * top(:,j)
! 1-2 break

    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_1d_2d_first!}}}

subroutine mul_op_1d_2d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(top, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  second * top(i,:)
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_1d_2d_second!}}}

subroutine mul_op_2d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % scalar

  do while (associated(temp_iter))
    second => second_iter % array

    ! do operation
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_2d_0d!}}}

subroutine mul_op_2d_1d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % array

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array
    jend = size(second, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  second(:,j) * top
! 1-2 break

    end do

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_2d_1d_first!}}}

subroutine mul_op_2d_1d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(second, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  second(i,:) * top
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_2d_1d_second!}}}

subroutine mul_op_2d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second * top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine mul_op_2d_2d!}}}



!***********************************************************************
! routine plus_operator
!
!> \brief Do plus on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> add them together, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine plus_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d
  op_name = '+'
! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were less than two operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer - 1) % number_of_dims
  b_nd = stack(stack_pointer) % number_of_dims

  ! call the right one
  if (a_nd == 0) then
    if (b_nd == 0) then

! 0d 1d
  call plus_op_0d_0d(stack, stack_pointer)
! 0d 1d

    else if (b_nd == 1) then

! 0d 2d
  call plus_op_0d_1d(stack, stack_pointer)
! 0d 2d

    else

! 1d 0d
  call plus_op_0d_2d(stack, stack_pointer)
! 1d 0d

    end if
  else if (a_nd == 1) then
    if (b_nd == 0) then

! 1d 1d same
  call plus_op_1d_0d(stack, stack_pointer)
! 1d 1d same

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 1d 1d diff
  call plus_op_1d_1d_same(stack, stack_pointer)
! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first
  call plus_op_1d_1d_diff(stack, stack_pointer)
! 1d 2d first

      end if
    else
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d2 % dimNames(1))) then

! 1d 2d second
  call plus_op_1d_2d_first(stack, stack_pointer)
! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d
  call plus_op_1d_2d_second(stack, stack_pointer)
! 2d 0d

      end if
    end if
  else 
    if (b_nd == 0) then

! 2d 1d first
  call plus_op_2d_0d(stack, stack_pointer)
! 2d 1d first

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d2 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 2d 1d second
  call plus_op_2d_1d_first(stack, stack_pointer)
! 2d 1d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
            trim(stack(stack_pointer) % d1 % dimNames(1))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 2d
  call plus_op_2d_1d_second(stack, stack_pointer)
! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end
  call plus_op_2d_2d(stack, stack_pointer)
! end

    end if
  end if
end subroutine plus_operator!}}}

subroutine plus_op_0d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d0, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d0

  ! do operation
  top => top_iter % scalar
  second => second_iter % scalar
  temp_iter % scalar = &

! 1-2 break
  second + top
! 1-2 break

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
end subroutine plus_op_0d_0d!}}}

subroutine plus_op_0d_1d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d0

  second => second_iter % scalar

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine plus_op_0d_1d!}}}

subroutine plus_op_0d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d0

  second => second_iter % scalar

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_0d_2d!}}}

subroutine plus_op_1d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d1

  top => top_iter % scalar 

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array 
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine plus_op_1d_0d!}}}

subroutine plus_op_1d_1d_same (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine plus_op_1d_1d_same!}}}

subroutine plus_op_1d_1d_diff (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, j, iend, jend

  ! allocate a temp for result
  call create_2d_field_from_1ds( &
    stack(stack_pointer - 1) % d1, stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array
  iend = size(second)

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top)

    do j = 1, jend
      do i = 1, iend
        temp_iter % array(i,j) = &

! 1-2 break
  second(i) + top(j)
! 1-2 break

      end do
    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_1d_1d_diff!}}}

subroutine plus_op_1d_2d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  second + top(:,j)
! 1-2 break

    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_1d_2d_first!}}}

subroutine plus_op_1d_2d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(top, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  second + top(i,:)
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_1d_2d_second!}}}

subroutine plus_op_2d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % scalar

  do while (associated(temp_iter))
    second => second_iter % array

    ! do operation
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_2d_0d!}}}

subroutine plus_op_2d_1d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % array

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array
    jend = size(second, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  second(:,j) + top
! 1-2 break

    end do

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_2d_1d_first!}}}

subroutine plus_op_2d_1d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(second, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  second(i,:) + top
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_2d_1d_second!}}}

subroutine plus_op_2d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second + top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine plus_op_2d_2d!}}}



!***********************************************************************
! routine minus_operator
!
!> \brief Do minus on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> subtract top from second, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine minus_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d
  op_name = '-'
! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were less than two operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer - 1) % number_of_dims
  b_nd = stack(stack_pointer) % number_of_dims

  ! call the right one
  if (a_nd == 0) then
    if (b_nd == 0) then

! 0d 1d
  call minus_op_0d_0d(stack, stack_pointer)
! 0d 1d

    else if (b_nd == 1) then

! 0d 2d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to subtract a 1d from a 0d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 0d 2d

    else

! 1d 0d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to subtract a 2d from a 0d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 1d 0d

    end if
  else if (a_nd == 1) then
    if (b_nd == 0) then

! 1d 1d same
  call minus_op_1d_0d(stack, stack_pointer)
! 1d 1d same

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 1d 1d diff
  call minus_op_1d_1d_same(stack, stack_pointer)
! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first
  call minus_op_1d_1d_diff(stack, stack_pointer)
! 1d 2d first

      end if
    else
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d2 % dimNames(1))) then

! 1d 2d second
  if (type_checking) then
    call mpas_log_write( &
      'Unable to subtract a 2d from a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to subtract a 2d from a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 2d 0d

      end if
    end if
  else 
    if (b_nd == 0) then

! 2d 1d first
  call minus_op_2d_0d(stack, stack_pointer)
! 2d 1d first

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d2 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 2d 1d second
  call minus_op_2d_1d_first(stack, stack_pointer)
! 2d 1d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
            trim(stack(stack_pointer) % d1 % dimNames(1))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 2d
  call minus_op_2d_1d_second(stack, stack_pointer)
! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end
  call minus_op_2d_2d(stack, stack_pointer)
! end

    end if
  end if
end subroutine minus_operator!}}}

subroutine minus_op_0d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d0, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d0

  ! do operation
  top => top_iter % scalar
  second => second_iter % scalar
  temp_iter % scalar = &

! 1-2 break
  second - top
! 1-2 break

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
end subroutine minus_op_0d_0d!}}}

subroutine minus_op_1d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d1

  top => top_iter % scalar 

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array 
    temp_iter % array = &

! 1-2 break
  second - top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine minus_op_1d_0d!}}}

subroutine minus_op_1d_1d_same (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second - top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine minus_op_1d_1d_same!}}}

subroutine minus_op_1d_1d_diff (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, j, iend, jend

  ! allocate a temp for result
  call create_2d_field_from_1ds( &
    stack(stack_pointer - 1) % d1, stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array
  iend = size(second)

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top)

    do j = 1, jend
      do i = 1, iend
        temp_iter % array(i,j) = &

! 1-2 break
  second(i) - top(j)
! 1-2 break

      end do
    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine minus_op_1d_1d_diff!}}}

subroutine minus_op_2d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % scalar

  do while (associated(temp_iter))
    second => second_iter % array

    ! do operation
    temp_iter % array = &

! 1-2 break
  second - top
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine minus_op_2d_0d!}}}

subroutine minus_op_2d_1d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % array

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array
    jend = size(second, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  second(:,j) - top
! 1-2 break

    end do

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine minus_op_2d_1d_first!}}}

subroutine minus_op_2d_1d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(second, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  second(i,:) - top
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine minus_op_2d_1d_second!}}}

subroutine minus_op_2d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  second - top
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine minus_op_2d_2d!}}}



!***********************************************************************
! routine div_operator
!
!> \brief Do div on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> divide the second by the top, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine div_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d
  op_name = '/'
! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were less than two operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer - 1) % number_of_dims
  b_nd = stack(stack_pointer) % number_of_dims

  ! call the right one
  if (a_nd == 0) then
    if (b_nd == 0) then

! 0d 1d
  call div_op_0d_0d(stack, stack_pointer)
! 0d 1d

    else if (b_nd == 1) then

! 0d 2d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to divide a 0d by a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 0d 2d

    else

! 1d 0d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to divide a 0d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 1d 0d

    end if
  else if (a_nd == 1) then
    if (b_nd == 0) then

! 1d 1d same
  call div_op_1d_0d(stack, stack_pointer)
! 1d 1d same

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 1d 1d diff
  call div_op_1d_1d_same(stack, stack_pointer)
! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first
  call div_op_1d_1d_diff(stack, stack_pointer)
! 1d 2d first

      end if
    else
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d2 % dimNames(1))) then

! 1d 2d second
  if (type_checking) then
    call mpas_log_write( &
      'Unable to divide a 1d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d
  if (type_checking) then
    call mpas_log_write( &
      'Unable to divide a 1d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM', MPAS_LOG_CRIT)
  end if
! 2d 0d

      end if
    end if
  else 
    if (b_nd == 0) then

! 2d 1d first
  call div_op_2d_0d(stack, stack_pointer)
! 2d 1d first

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d2 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 2d 1d second
  call div_op_2d_1d_first(stack, stack_pointer)
! 2d 1d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
            trim(stack(stack_pointer) % d1 % dimNames(1))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 2d
  call div_op_2d_1d_second(stack, stack_pointer)
! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end
  call div_op_2d_2d(stack, stack_pointer)
! end

    end if
  end if
end subroutine div_operator!}}}

function safe_divide_0d_0d(second, top)
  implicit none
  real (kind=RKIND), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND) :: safe_divide_0d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_0d_0d = second / top
  else
    safe_divide_0d_0d = huge(second)
  end if
end function safe_divide_0d_0d

function safe_divide_1d_0d(second, top)
  implicit none
  real (kind=RKIND), dimension(:), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND), dimension(size(second)) :: safe_divide_1d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_1d_0d = second / top
  else
    safe_divide_1d_0d = huge(second)
  end if
end function safe_divide_1d_0d

function safe_divide_2d_0d(second, top)
  implicit none
  real (kind=RKIND), dimension(:, :), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND), dimension(size(second, 1), size(second, 2)) :: &
    safe_divide_2d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_2d_0d = second / top
  else
    safe_divide_2d_0d = huge(second)
  end if
end function safe_divide_2d_0d

function safe_divide_1d_1d(second, top)
  implicit none
  real (kind=RKIND), dimension(:), intent(in) :: second
  real (kind=RKIND), dimension(:), intent(in) :: top

  real (kind=RKIND), dimension(size(second)) :: safe_divide_1d_1d

  where (abs(top) > 0.0_RKIND)
    safe_divide_1d_1d = second / top
  elsewhere
    safe_divide_1d_1d = huge(second)
  end where
end function safe_divide_1d_1d

function safe_divide_2d_2d(second, top)
  implicit none
  real (kind=RKIND), dimension(:, :), intent(in) :: second
  real (kind=RKIND), dimension(:, :), intent(in) :: top

  real (kind=RKIND), dimension(size(second, 1), size(second, 2)) :: &
    safe_divide_2d_2d

  where (abs(top) > 0.0_RKIND)
    safe_divide_2d_2d = second / top
  elsewhere
    safe_divide_2d_2d = huge(second)
  end where
end function safe_divide_2d_2d

subroutine div_op_0d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d0, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d0

  ! do operation
  top => top_iter % scalar
  second => second_iter % scalar
  temp_iter % scalar = &

! 1-2 break
  safe_divide_0d_0d(second, top)
! 1-2 break

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
end subroutine div_op_0d_0d!}}}

subroutine div_op_1d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d1

  top => top_iter % scalar 

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array 
    temp_iter % array = &

! 1-2 break
  safe_divide_1d_0d(second, top)
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine div_op_1d_0d!}}}

subroutine div_op_1d_1d_same (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  safe_divide_1d_1d(second, top)
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine div_op_1d_1d_same!}}}

subroutine div_op_1d_1d_diff (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, j, iend, jend

  ! allocate a temp for result
  call create_2d_field_from_1ds( &
    stack(stack_pointer - 1) % d1, stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array
  iend = size(second)

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top)

    do j = 1, jend
      do i = 1, iend
        temp_iter % array(i,j) = &

! 1-2 break
  safe_divide_0d_0d(second(i), top(j))
! 1-2 break

      end do
    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine div_op_1d_1d_diff!}}}

subroutine div_op_2d_0d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % scalar

  do while (associated(temp_iter))
    second => second_iter % array

    ! do operation
    temp_iter % array = &

! 1-2 break
  safe_divide_2d_0d(second, top)
! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine div_op_2d_0d!}}}

subroutine div_op_2d_1d_first (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  top => top_iter % array

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array
    jend = size(second, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
  safe_divide_1d_1d(second(:,j), top)
! 1-2 break

    end do

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine div_op_2d_1d_first!}}}

subroutine div_op_2d_1d_second (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second
  integer :: i, iend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array

    iend = size(second, 1)

    do i = 1, iend
      temp_iter % array(i,:) = &

! 1-2 break
  safe_divide_1d_1d(second(i,:), top)
! 1-2 break

    end do

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine div_op_2d_1d_second!}}}

subroutine div_op_2d_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field2DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:,:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d2

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    second => second_iter % array
    temp_iter % array = &

! 1-2 break
  safe_divide_2d_2d(second, top)
! 1-2 break

    top_iter => top_iter % next
    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
end subroutine div_op_2d_2d!}}}

!***********************************************************************
! routine sum_operator
!
!> \brief Do sum on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take sum argument off the stack and
!> sum along the first dimension, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine sum_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd
  character (len=StrKIND) ::  op_name

! start -> 1d
  op_name = 'sum'
! start -> 1d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 1) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were no operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer) % number_of_dims

  if (type_checking) then
    if (a_nd < 1) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when the ' // &
        'operand on the stack is 0d, in the RPN calculator AM')
    end if
  end if

  ! call the right one
  if (a_nd == 1) then

! 1d -> 2d
  call sum_op_1d(stack, stack_pointer)
! 1d -> 2d

  else

! 2d -> end
  call sum_op_2d(stack, stack_pointer)
! 2d -> end

  end if
end subroutine sum_operator

subroutine sum_op_1d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp
  type (field1DReal), pointer :: top_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), pointer :: reduced

  ! allocate a temp for result
  call create_0d_field_from_1d(stack(stack_pointer) % d1, temp)

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1

  ! initial value
  reduced => temp % scalar
  temp % scalar = &

! 1-2 break
  0
! 1-2 break

  do while (associated(top_iter))
    ! do operation
    top => top_iter % array

    temp % scalar = &

! 2-3 break
  reduced + sum(top)
! 2-3 break

    top_iter => top_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  ! set stack
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
end subroutine sum_op_1d!}}}

subroutine sum_op_2d (stack, stack_pointer)!{{{
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: reduced
  integer :: j, jend

  ! allocate a temp for result
  call create_1d_field_from_2d(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2

  ! initial value
  reduced => temp_iter % array
  temp_iter % array = &

! 1-2 break
  0
! 1-2 break

  do while (associated(top_iter))
    ! do operation
    top => top_iter % array
    reduced => temp_iter % array

    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(j) = &

! 2-3 break
  reduced(j) + sum(top(:,j))
! 2-3 break

    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  ! set stack
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
end subroutine sum_op_2d!}}}

end module ocn_rpn_calculator
! vim: foldmethod=marker
