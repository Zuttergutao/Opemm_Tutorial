from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout



# 输入文件
profile="protein.pdb"
ligfile="lig.pdb"
ligtop="lig.xml"
ph=7.0

# 系统设置
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
constraints = HBonds
rigidWater = True
hydrogenMass = 1.5*amu
ewaldErrorTolerance = 0.0005
constraintTolerance = 0.000001

# 积分选项
dt = 0.002*picoseconds
temperature = 306*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 100

# 模拟选项
steps = 100000
equilibrationSteps = 1000
platform = Platform.getPlatformByName('CUDA')
recordinpdb=False
recordindcd=True
stdreport=True
pdbReporter = PDBReporter("trajectory.pdb",5000)
dcdReporter = DCDReporter('trajectory.dcd', 5000)
stdReporter = StateDataReporter(stdout, 1000, totalSteps=steps,
    step=True, speed=True, progress=True,kineticEnergy=True, potentialEnergy=True, temperature=True, separator='\t')
dataReporter = StateDataReporter("log.txt", 1000, totalSteps=steps,
    step=True, speed=True, progress=True,kineticEnergy=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('checkpoint.chk', 50000)

# 处理蛋白
print("处理蛋白中--------------------------------")
fixer=PDBFixer(filename=profile)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(True)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(ph)
PDBFile.writeFile(fixer.topology,fixer.positions,open("protein_fixed.pdb","w"))
print("蛋白处理完成------------------------------")

# 合并蛋白和配体并生成盒子
print("合并蛋白-配体并生成模拟体系---------------")
lig=PDBFile(ligfile)
pro=PDBFile("protein_fixed.pdb")
modeller=Modeller(pro.topology,pro.positions)
modeller.add(lig.topology,lig.positions)
forcefield=ForceField("amber14-all.xml","amber14/tip3p.xml",ligtop)
modeller.addSolvent(padding=8*angstrom,forcefield=forcefield,neutralize=True)
mergedTopology=modeller.topology
mergedPositions=modeller.positions
PDBFile.writeFile(mergedTopology,mergedPositions,open("complex.pdb","w"))


# 构建模拟系统
print("模拟准备中--------------------------------")
system=forcefield.createSystem(mergedTopology,nonbondedMethod=nonbondedMethod,nonbondedCutoff=nonbondedCutoff,
        constraints=constraints,rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance,hydrogenMass=hydrogenMass)
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(mergedTopology, system, integrator, platform)
simulation.context.setPositions(mergedPositions)


# 最小化
print("能量最小化--------------------------------")
simulation.minimizeEnergy()

# 预平衡
print("平衡中------------------------------------")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# 成品模拟
print("Simulating--------------------------------")
if recordinpdb:
    simulation.reporters.append(pdbReporter)
if recordindcd:
    simulation.reporters.append(dcdReporter)
if stdreport:
    simulation.reporters.append(stdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)
