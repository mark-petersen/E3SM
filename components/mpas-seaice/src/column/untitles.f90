real (kind=dbl_kind), intent(in) :: &
initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

real (kind=dbl_kind), intent(in) :: &
ratio_Si2N_diatoms, &   ! algal Si to N (mol/mol)
ratio_Si2N_sp     , &
ratio_Si2N_phaeo  , &
ratio_S2N_diatoms , &   ! algal S  to N (mol/mol)
ratio_S2N_sp      , &
ratio_S2N_phaeo   , &
ratio_Fe2C_diatoms, &   ! algal Fe to C  (umol/mol)
ratio_Fe2C_sp     , &
ratio_Fe2C_phaeo  , &
ratio_Fe2N_diatoms, &   ! algal Fe to N  (umol/mol)
ratio_Fe2N_sp     , &
ratio_Fe2N_phaeo  , &
ratio_Fe2DON      , &   ! Fe to N of DON (nmol/umol)
ratio_Fe2DOC_s    , &   ! Fe to C of DOC (nmol/umol) saccharids
ratio_Fe2DOC_l    , &   ! Fe to C of DOC (nmol/umol) lipids 
tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
chlabs_diatoms   , & ! chl absorption (1/m/(mg/m^3))
chlabs_sp        , & !
chlabs_phaeo     , & !
alpha2max_low_diatoms , & ! light limitation (1/(W/m^2))  
alpha2max_low_sp      , & 
alpha2max_low_phaeo   , & 
beta2max_diatoms , & ! light inhibition (1/(W/m^2))  
beta2max_sp      , & 
beta2max_phaeo   , & 
mu_max_diatoms   , & ! maximum growth rate (1/day)       
mu_max_sp        , & 
mu_max_phaeo     , & 
grow_Tdep_diatoms, & ! Temperature dependence of growth (1/C)
grow_Tdep_sp     , & 
grow_Tdep_phaeo  , & 
fr_graze_diatoms , & ! Fraction grazed
fr_graze_sp      , & 
fr_graze_phaeo   , & 
mort_pre_diatoms , & ! Mortality (1/day)
mort_pre_sp      , & 
mort_pre_phaeo   , & 
mort_Tdep_diatoms, & ! T dependence of mortality (1/C)
mort_Tdep_sp     , &  
mort_Tdep_phaeo  , &  
k_exude_diatoms  , & ! algal exudation (1/d)
k_exude_sp       , &  
k_exude_phaeo    , &  
K_Nit_diatoms    , & ! nitrate half saturation (mmol/m^3)
K_Nit_sp        , &  
K_Nit_phaeo      , &  
K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
K_Am_sp         , &   
K_Am_phaeo       , &   
K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
K_Sil_sp        , &   
K_Sil_phaeo      , &   
K_Fe_diatoms     , & ! iron half saturation (nM)
K_Fe_sp         , &   
K_Fe_phaeo       , &    
f_don_protein    , & ! fraction of spilled grazing to proteins          
kn_bac_protein   , & ! Bacterial degredation of DON (1/d)               
f_don_Am_protein , & ! fraction of remineralized DON to ammonium        
f_doc_s         , & ! fraction of mortality to DOC 
f_doc_l         , &   
f_exude_s        , & ! fraction of exudation to DOC
f_exude_l        , & 
k_bac_s         , & ! Bacterial degredation of DOC (1/d)
k_bac_l         , & 
algaltype_diatoms  , & ! mobility type
algaltype_sp       , & !
algaltype_phaeo    , & !
nitratetype        , & !
ammoniumtype       , & !
silicatetype       , & !
dmspptype         , & !
dmspdtype         , & !
humtype           , & !
doctype_s         , & !
doctype_l         , & !
dictype_1         , & !
dontype_protein    , & !
fedtype_1         , & !
feptype_1         , & !
zaerotype_bc1      , & !
zaerotype_bc2      , & !
zaerotype_dust1    , & !
zaerotype_dust2    , & !
zaerotype_dust3    , & !
zaerotype_dust4    , & !
ratio_C2N_diatoms  , & ! algal C to N ratio (mol/mol)
ratio_C2N_sp       , & !
ratio_C2N_phaeo    , & !
ratio_chl2N_diatoms, & ! algal chlorophyll to N ratio (mg/mmol)
ratio_chl2N_sp     , & !
ratio_chl2N_phaeo  , & !
F_abs_chl_diatoms  , & ! scales absorbed radiation for dEdd
F_abs_chl_sp       , & !
F_abs_chl_phaeo    , & !
ratio_C2N_proteins     ! ratio of C to N in proteins (mol/mol)   
