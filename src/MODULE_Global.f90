
!************************************************************
!  File Type:     MODULE
!  File Content:  Modules containing global variables and setup
!
!************************************************************


   MODULE GlobVar                                  ! Global variables used MANY places
      Use dpDef

      Real(dp) :: Ptot,Temp                        ! Ptot is pressure in bara and Temp is temperature in Kelvin.
      Real(dp), Parameter :: R     =8.314473_dp    ! Gas constant for thermodynamic calculations. Units = [J/K mole]
      Real(dp), parameter :: Rgas = 8.314473_dp

   End MODULE GlobVar




   MODULE Par

      ! ManualUpdate check that these components are updated

      Integer, Parameter :: NoAq=13          ! Number of species in aqueous model
      Integer, Parameter :: NoComp=20        ! Number of defined components in CO2 - gas components and possible solids
      Integer, Parameter :: NoHC=20          ! Gas components - included in EOS
      Integer, Parameter :: Nsolid=1         ! Number of solids
      Integer, Parameter :: NInputComp=20    ! Number of input components

      Integer NGasComp, NSolComp

  !    Integer, Parameter :: MaxNsalt=22      ! Number of salts included in model
  !    Integer, Parameter :: MaxLev=2         ! Number of levels in calculation 1=Before precipiation, 2=Equilibrium
  !    Integer, Parameter :: MaxPoints=100    ! Maximum number of points in a profile
  !    Integer, Parameter :: NoPseudo=13       ! Number of pseudocomponents (maksiumum)


   End MODULE Par





   MODULE ResHandling
      Integer :: UserChoices(0:10)

   End MODULE ResHandling


   MODULE Path

      Character Paths(4)*500

      Character FileName*500

      Integer, Parameter :: InFile=19
      Integer, Parameter :: OutFile=20
      Integer, Parameter :: LogFilePath=21
      Integer, Parameter :: LogFileCalc=22
      Integer, Parameter :: FlashResFile=23

      Integer, Parameter :: DataPath      = 1
      Integer, Parameter :: SettingsPath  = 2
      Integer, Parameter :: TempPath      = 3

      Logical LOG

   End MODULE Path



   MODULE Files

      Logical ReadFromFile   ! Determines if all the data should be read from file or set directly

      Character, Parameter :: AqSetupFile*13          = 'AqSpecies.dat'   ! Aqeous species, charge, Mw etc.
      !Character, Parameter :: HCSetupFile*13          = 'HCSpecies.dat'
      Character, Parameter :: EqSEtupFile*14          = 'Equilibria.dat'
      !Character, Parameter :: SimulationInputFile*19  = 'SimulationInput.dat'
      Character, Parameter :: EqThermFile*12          = 'Eq_therm.dat'
      Character, Parameter :: WaterDataFile*13        = 'Waterdata.dat'
      !Character, Parameter :: SaltDataFile*12         = 'Saltdata.dat'
      !Character, Parameter :: LogFilePathFileName*15  = 'LogFilePath.dat'
      !Character, Parameter :: LogFileCalcFileName*15  = 'LogFileCalc.dat'
      Character, Parameter :: VolFileName*8           = 'Vols.dat'
      Character, Parameter :: EOS_PureComp*11         = 'EOSData.dat'
      Character, Parameter :: EOS_BIN_File*10         = 'EOSBin.dat'

      Character, Parameter :: GasSolubFile*12         = 'GasSolub.dat'
      !Character, Parameter :: BinParFile*10          = 'Binpar.dat'
      !Character, Parameter :: FlashRes*12            = 'FlashRes.dat'
      Character, Parameter :: GasSetupFile*14         = 'GasSpecies.dat'
      Character, Parameter :: CompThermoFile*14       ='CompThermo.dat'
      Character, Parameter :: ReactionsFile*13        ='Reactions.dat'

   End MODULE Files


   MODULE ErrMSGS
      Integer Info
      Character ErrMSG*500

   End MODULE ErrMSGS