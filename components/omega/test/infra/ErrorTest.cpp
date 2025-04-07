//===-- Test driver for OMEGA error handler ----------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA error handler
///
/// This driver tests the error handling capabilities for the OMEGA
/// model. In particular, it tests error checking for various error types
/// and the writing of stack traces for critical errors.
///
//
//===-----------------------------------------------------------------------===/

#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Create a call tree for less critical errors to test passing return codes
// up the chain

// At lowest level, create and return a warning message
Error returnErrorRoutine3(void) {
   Error RetVal; // init to success
   int Level = 3;
   ERROR_RETURN(RetVal, ErrorCode::Warn, "Test warning at level {}", Level);
}

// At second level, add another warning message
Error returnErrorRoutine2(void) {
   Error RetVal; // init to success
   int Level = 2;
   RetVal += returnErrorRoutine3(); // tests accumulator
   ERROR_RETURN(RetVal, ErrorCode::Warn, "Test warning at level {}", Level);
}

Error returnErrorRoutine1(void) {
   Error ErrorTmp1; // init to success
   int Level       = 1;
   Error ErrorTmp2 = returnErrorRoutine2();
   Error RetVal    = ErrorTmp1 + ErrorTmp2;
   // upgrade from warning to fail
   ERROR_RETURN(RetVal, ErrorCode::Fail, "Test error at level {}", Level);
}

//------------------------------------------------------------------------------
// Create critical error call chain to test an abort and a correct stack trace

void passError2(Error &Err) {
   R8 MsgArgReal         = 1.23456789012;
   I4 MsgArgInt          = 5;
   std::string MsgArgStr = "Random";
   if (!Err.isSuccess())
      ERROR_CRITICAL("Expected critical error with args: {} {} {}", MsgArgStr,
                     MsgArgInt, MsgArgReal);
   // Alternative critical failure tests
   // ERROR_CHECK_CRITICAL(Err, "Expected critical error with args: {} {} {}",
   //                     MsgArgStr, MsgArgInt, MsgArgReal);
   // ERROR_REQUIRE(!Err.isSuccess(),
   //              "Expected Error with args: {} {} {}", MsgArgStr, MsgArgInt,
   //               MsgArgReal);
   // This one requires a debug build
   // ERROR_ASSERT(!Err.isSuccess(),
   //             "Expected Error with args: {} {} {}", MsgArgStr, MsgArgInt,
   //             MsgArgReal);
}

void passError1(Error &Err) { passError2(Err); }

//------------------------------------------------------------------------------
// Main test driver
int main(int argc, char **argv) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   initLogging(DefEnv);

   LOG_INFO("------ Error handler test ------");
   // Test default constructor and success check
   Error ErrDefault;
   if (ErrDefault.isSuccess()) {
      LOG_INFO("Default error constructor: PASS");
   } else {
      LOG_ERROR("Default error constructor: FAIL");
   }

   // Test constructors with no message and boolean tests
   Error ErrorTestSuccess2(ErrorCode::Success);
   if (ErrorTestSuccess2.isSuccess() and !ErrorTestSuccess2.isFail()) {
      LOG_INFO("No message error constructor success: PASS");
   } else {
      LOG_ERROR("No message error constructor success: FAIL");
   }
   Error ErrorTestWarn2(ErrorCode::Warn);
   if (!ErrorTestWarn2.isSuccess() and !ErrorTestWarn2.isFail() and
       ErrorTestWarn2.Code == ErrorCode::Warn) {
      LOG_INFO("No message error constructor warn: PASS");
   } else {
      LOG_ERROR("No message error constructor warn: FAIL");
   }
   Error ErrorTestFail2(ErrorCode::Fail);
   if (ErrorTestFail2.isFail() and !ErrorTestFail2.isSuccess()) {
      LOG_INFO("No message error constructor fail: PASS");
   } else {
      LOG_ERROR("No message error constructor fail: FAIL");
   }

   // Test constructor with message
   // Create a total error to accumulate codes and messages
   std::string ExpectMsg =
       "   [error] [ErrorTest.cpp:121] Initialized to fail String 1\n";
   Error TotalError(ErrorCode::Fail, __LINE__, __FILE__,
                    "Initialized to fail {} {}", "String", 1);
   if (TotalError.isFail() and TotalError.Msg == ExpectMsg) {
      LOG_INFO("Error constructor with message: PASS");
   } else {
      LOG_ERROR("Error constructor with message: FAIL");
      LOG_ERROR(" Error constructor with message expected: {}", ExpectMsg);
      LOG_ERROR(" Error constructor with message actual:   {}", TotalError.Msg);
   }

   // Test reset of error code
   TotalError.reset();
   if (TotalError.isSuccess()) {
      LOG_INFO("Reset of error code: PASS");
   } else {
      ERROR_CHECK(TotalError, "Reset of error code: FAIL");
   }

   // Call routines to check error return functions
   TotalError += returnErrorRoutine1();

   // Test ERROR_CHECK for adding and printing error message
   ERROR_CHECK(TotalError, "Encountered error in returnErrorRoutine1");

   // Test ERROR_CHECK for adding warning and error message
   ERROR_CHECK_WARN(TotalError, "Encountered error in returnErrorRoutine1");

   // The above calls have reset the error so create a new failure
   TotalError += Error(ErrorCode::Fail, __LINE__, __FILE__, "Error in driver");

   // Test critical error
   passError1(TotalError);

   // Finalize environments
   MPI_Finalize();

   return RetVal;
}
//------------------------------------------------------------------------------
