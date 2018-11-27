SetCurrentWorkDirectory (FilePath = u'C:\\Users\\cma86\\AppData\\Local\\Temp\\WB_CBMOONSHOT1C16_cma86_6232_2\\unsaved_project_files\\dp0\\Study\\AIM\\physics2', isDifferentContents = True)
Root.Analysis['DefaultAnalysis'].SetState (State = {'AnalysisType': {'Option': 'Steady'}})
Root.Analysis['DefaultAnalysis'].SetState (State = {'SolutionControls': {}})
Root.Analysis['DefaultAnalysis'].SetState (State = {'SolverName': 'FLUX'})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].SetState (State = {'PhysicsRegionLocation': []})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].SetState (State = {'PhysicsType': {'PhysicsRegionType': 'Fluid', 'PhysicsTypeList': ['Fluid']}})
SynchroniseTopology (State = '2\n0\n0\n1\n1 0 40\n10\n11 0 189\n1 0 76\n2 0 77\n3 0 78\n4 0 79\n5 0 80\n6 0 81\n7 0 82\n8 0 83\n10 0 188\n24\n22 0 190\n24 0 192\n21 0 75\n26 0 194\n1 0 55\n2 0 56\n25 0 193\n4 0 58\n5 0 59\n6 0 60\n7 0 61\n8 0 62\n9 0 63\n10 0 64\n11 0 65\n12 0 66\n13 0 67\n14 0 68\n15 0 69\n16 0 70\n17 0 71\n18 0 72\n23 0 191\n20 0 74\n16\n15 0 182\n16 0 183\n18 0 187\n20 0 184\n1 0 41\n2 0 42\n19 0 185\n5 0 45\n6 0 46\n7 0 47\n8 0 48\n9 0 49\n10 0 50\n11 0 51\n17 0 186\n13 0 53\n106\n32 1 11 0\n32 1 1 0\n32 1 2 0\n32 1 3 0\n32 1 4 0\n32 1 5 0\n32 1 6 0\n32 1 7 0\n32 1 8 0\n32 1 10 0\n21 11 22 1\n21 11 24 0\n21 11 21 0\n21 11 26 1\n10 22 15 0\n10 22 16 1\n10 24 15 0\n10 24 18 1\n10 21 18 0\n10 21 20 1\n10 26 16 0\n10 26 20 1\n21 1 1 0\n21 1 2 1\n21 1 25 0\n21 1 26 0\n21 1 4 0\n10 1 1 0\n10 1 2 1\n10 2 19 0\n10 2 2 1\n10 25 19 0\n10 25 16 1\n10 4 20 0\n10 4 1 1\n21 2 5 0\n21 2 6 1\n21 2 7 1\n21 2 2 0\n10 5 2 0\n10 5 5 1\n10 6 6 0\n10 6 5 1\n10 7 19 0\n10 7 6 1\n21 3 8 0\n21 3 9 1\n21 3 10 0\n21 3 6 0\n10 8 5 0\n10 8 7 1\n10 9 8 0\n10 9 7 1\n10 10 8 0\n10 10 6 1\n21 4 11 0\n21 4 12 1\n21 4 13 1\n21 4 9 0\n10 11 7 0\n10 11 9 1\n10 12 10 0\n10 12 9 1\n10 13 8 0\n10 13 10 1\n21 5 14 0\n21 5 15 1\n21 5 16 0\n21 5 12 0\n10 14 9 0\n10 14 11 1\n10 15 17 0\n10 15 11 1\n10 16 17 0\n10 16 10 1\n21 6 17 0\n21 6 18 1\n21 6 24 1\n21 6 23 1\n21 6 15 0\n10 17 11 0\n10 17 13 1\n10 18 18 0\n10 18 13 1\n10 23 17 0\n10 23 15 1\n21 7 20 0\n21 7 4 1\n21 7 21 1\n21 7 18 0\n10 20 13 0\n10 20 1 1\n21 8 1 1\n21 8 20 1\n21 8 17 1\n21 8 14 1\n21 8 11 1\n21 8 8 1\n21 8 5 1\n21 10 22 0\n21 10 25 1\n21 10 7 0\n21 10 10 1\n21 10 13 0\n21 10 16 1\n21 10 23 0\n')
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].SetState (State = {'PhysicsType': {'PhysicsRegionType': 'Fluid', 'PhysicsTypeList': ['Fluid']}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].SetState (State = {'PhysicsRegionLocation': ['BODY1']})
Root.Library.MaterialsDB.MaterialBehavior['MaterialAtState 1'].SetState (State = {'MaterialModel': 'ConstantProperties',
 'MaterialName': 'Material 1',
 'PhysicalState': 'Liquid'})
Root.Library.MaterialsDB.Material['Material 1'].SetState (State = {'LiquidProperties': {'Density': {'Option': 'Value',
                                  'PropertyValue': '997.04763676 [kg m^-3]'},
                      'DynamicViscosity': {'Option': 'Value',
                                           'PropertyValue': '0.000890022489078 [Pa s]'},
                      'SpecificHeatCapacityCp': {'Option': 'Value',
                                                 'PropertyValue': '4181.31499079 [J kg^-1 K^-1]'},
                      'ThermalConductivity': {'Option': 'Value',
                                              'PropertyValue': '0.60651608022 [W m^-1 K^-1]'}},
 'MolarMass': '0.018015268 [kg mol^-1]'})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].SetState (State = {'InitialConditions': {},
 'NumericalControls': {},
 'PhaseDefinitions': {},
 'PhaseModels': {},
 'PhysicsRegionReferenceFrame': {},
 'VariableControls': {}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].PhaseDefinitions.PhaseDefinition['MaterialAssignment 1'].SetState (State = {'MaterialBehaviorName': 'MaterialAtState 1',
 'MaterialDistribution': 'Continuous'})
Root.Library.MaterialsDB.MaterialBehavior['MaterialAtState 1'].SetState (State = {'MaterialModel': 'ConstantProperties',
 'MaterialName': 'Material 1',
 'PhysicalState': 'Liquid'})
Root.Library.MaterialsDB.Material['Material 1'].SetState (State = {'LiquidProperties': {'Density': {'Option': 'Value',
                                  'PropertyValue': '997.04763676 [kg m^-3]'},
                      'DynamicViscosity': {'Option': 'Value',
                                           'PropertyValue': '0.000890022489078 [Pa s]'},
                      'SpecificHeatCapacityCp': {'Option': 'Value',
                                                 'PropertyValue': '4181.31499079 [J kg^-1 K^-1]'},
                      'ThermalConductivity': {'Option': 'Value',
                                              'PropertyValue': '0.60651608022 [W m^-1 K^-1]'}},
 'MolarMass': '0.018015268 [kg mol^-1]'})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].PhaseDefinitions.PhaseDefinition['MaterialAssignment 1'].SetState (State = {'MaterialBehaviorName': 'MaterialAtState 1',
 'MaterialDistribution': 'Continuous'})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 1'].SetState (State = {})
SetState (ObjectPath = '/FLUX', State = {'Analysis:DefaultAnalysis': {'PhysicsRegion:PhysicsRegion 2': {'Boundary:Boundary 1': {'BoundaryType': 'Inlet'}}}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 1'].SetState (State = {'BoundaryLocations': ['FACE7']})
SetState (ObjectPath = '/FLUX', State = {'Analysis:DefaultAnalysis': {'PhysicsRegion:PhysicsRegion 2': {'Boundary:Boundary 1': {'Flow': {'Velocity': {'Magnitude': '0.01 [m s^-1]'}}}}}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 2'].SetState (State = {})
SetState (ObjectPath = '/FLUX', State = {'Analysis:DefaultAnalysis': {'PhysicsRegion:PhysicsRegion 2': {'Boundary:Boundary 2': {'BoundaryType': 'Outlet'}}}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 2'].SetState (State = {'BoundaryLocations': ['FACE4']})
SetState (ObjectPath = '/FLUX', State = {'Analysis:DefaultAnalysis': {'PhysicsRegion:PhysicsRegion 2': {'Boundary:Boundary 2': {'Flow': {'Pressure': {'GaugeStaticPressure': '0 [Pa]'}}}}}})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 3'].SetState (State = {})
Root.Analysis['DefaultAnalysis'].PhysicsRegion['PhysicsRegion 2'].Boundary['Boundary 3'].SetState (State = {'BoundaryLocations': ['FACE1',
                       'FACE2',
                       'FACE3',
                       'FACE5',
                       'FACE6',
                       'FACE8',
                       'FACE10',
                       'FACE11']})
SetState (ObjectPath = '/FLUX', State = {'Analysis:DefaultAnalysis': {'PhysicsRegion:PhysicsRegion 2': {'Boundary:Boundary 3': {'BoundaryType': 'Wall'}}}})
SynchroniseTopology (State = '2\n0\n0\n1\n1 0 40\n10\n11 0 189\n1 0 76\n2 0 77\n3 0 78\n4 0 79\n5 0 80\n6 0 81\n7 0 82\n8 0 83\n10 0 188\n24\n22 0 190\n24 0 192\n21 0 75\n26 0 194\n1 0 55\n2 0 56\n25 0 193\n4 0 58\n5 0 59\n6 0 60\n7 0 61\n8 0 62\n9 0 63\n10 0 64\n11 0 65\n12 0 66\n13 0 67\n14 0 68\n15 0 69\n16 0 70\n17 0 71\n18 0 72\n23 0 191\n20 0 74\n16\n15 0 182\n16 0 183\n18 0 187\n20 0 184\n1 0 41\n2 0 42\n19 0 185\n5 0 45\n6 0 46\n7 0 47\n8 0 48\n9 0 49\n10 0 50\n11 0 51\n17 0 186\n13 0 53\n106\n32 1 11 0\n32 1 1 0\n32 1 2 0\n32 1 3 0\n32 1 4 0\n32 1 5 0\n32 1 6 0\n32 1 7 0\n32 1 8 0\n32 1 10 0\n21 11 22 1\n21 11 24 0\n21 11 21 0\n21 11 26 1\n10 22 15 0\n10 22 16 1\n10 24 15 0\n10 24 18 1\n10 21 18 0\n10 21 20 1\n10 26 16 0\n10 26 20 1\n21 1 1 0\n21 1 2 1\n21 1 25 0\n21 1 26 0\n21 1 4 0\n10 1 1 0\n10 1 2 1\n10 2 19 0\n10 2 2 1\n10 25 19 0\n10 25 16 1\n10 4 20 0\n10 4 1 1\n21 2 5 0\n21 2 6 1\n21 2 7 1\n21 2 2 0\n10 5 2 0\n10 5 5 1\n10 6 6 0\n10 6 5 1\n10 7 19 0\n10 7 6 1\n21 3 8 0\n21 3 9 1\n21 3 10 0\n21 3 6 0\n10 8 5 0\n10 8 7 1\n10 9 8 0\n10 9 7 1\n10 10 8 0\n10 10 6 1\n21 4 11 0\n21 4 12 1\n21 4 13 1\n21 4 9 0\n10 11 7 0\n10 11 9 1\n10 12 10 0\n10 12 9 1\n10 13 8 0\n10 13 10 1\n21 5 14 0\n21 5 15 1\n21 5 16 0\n21 5 12 0\n10 14 9 0\n10 14 11 1\n10 15 17 0\n10 15 11 1\n10 16 17 0\n10 16 10 1\n21 6 17 0\n21 6 18 1\n21 6 24 1\n21 6 23 1\n21 6 15 0\n10 17 11 0\n10 17 13 1\n10 18 18 0\n10 18 13 1\n10 23 17 0\n10 23 15 1\n21 7 20 0\n21 7 4 1\n21 7 21 1\n21 7 18 0\n10 20 13 0\n10 20 1 1\n21 8 1 1\n21 8 20 1\n21 8 17 1\n21 8 14 1\n21 8 11 1\n21 8 8 1\n21 8 5 1\n21 10 22 0\n21 10 25 1\n21 10 7 0\n21 10 10 1\n21 10 13 0\n21 10 16 1\n21 10 23 0\n')
UpdateMesh (Port = 62466, Timeout = 1)
InitializeSolution ()
Solve ()
CreateAllResultDomains ()
CreateWireframe (Name = 'DefaultWireframe')
CreateSurfaceUniformSampledPoints (Name = 'Surface seed points for Result 70 on FACE8', LocationList = ['Topology:FACE8'], NumberOfSamples = 20)
CreateVolumeStreamLine (Name = 'StreamLine 1', SeedList = ['SampledPoints:Surface seed points for Result 70 on FACE8'], MaxSteps = 1000, Direction = 'Forward', SkipFactor = 1, StepSizeControl = 'Automatic')
CreateResult (Name = 'Result PathlineTime for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'PathlineTime')
CreateResult (Name = 'Result StreamlineVelocity for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'Velocity')
CreateResult (Name = 'Result StreamlineVorticity for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'Vorticity')
CreateSurfaceUniformSampledPoints (Name = 'Surface seed points for Result 70 on FACE8', LocationList = ['Topology:FACE8'], NumberOfSamples = 100)
CreateVolumeStreamLine (Name = 'StreamLine 1', SeedList = ['SampledPoints:Surface seed points for Result 70 on FACE8'], MaxSteps = 1000, Direction = 'Forward', SkipFactor = 1, StepSizeControl = 'Automatic')
CreateResult (Name = 'Result PathlineTime for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'PathlineTime')
CreateResult (Name = 'Result StreamlineVelocity for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'Velocity')
CreateResult (Name = 'Result StreamlineVorticity for StreamLine 1', Location = 'Streamline:StreamLine 1', Variable = 'Vorticity')
CreateVolumeGridSampledPoints (Name = 'Volume seed points for Result 71 on BODY1', NumberOfSamples = 100, LocationList = ['Topology:BODY1'])
CreateResultsForVariable (NameList = ['Result 71 for Volume seed points for Result 71 on BODY1'], LocationList = ['VolumeSampledPoints:Volume seed points for Result 71 on BODY1'], Variable = 'Velocity')
DeletePostObject (ObjectPath = 'Result:Result 71 for Volume seed points for Result 71 on BODY1')
DeletePostObject (ObjectPath = 'VolumeSampledPoints:Volume seed points for Result 71 on BODY1')
CreateSurfaceUniformSampledPoints (Name = 'Surface seed points for Result 71 on Boundary 2', LocationList = ['Boundary:Boundary 2'], NumberOfSamples = 100)
CreateResultsForVariable (NameList = ['Result 71 for Surface seed points for Result 71 on Boundary 2'], LocationList = ['SampledPoints:Surface seed points for Result 71 on Boundary 2'], Variable = 'Velocity')
DeletePostObject (ObjectPath = 'Result:Result 71 for Surface seed points for Result 71 on Boundary 2')
DeletePostObject (ObjectPath = 'SampledPoints:Surface seed points for Result 71 on Boundary 2')
