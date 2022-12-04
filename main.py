# -*- coding: mbcs -*-
from part import *
import regionToolset
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import abaqus as aq

import os
from os.path import join
import json

    
directory = "C:\\Users\\Nikhil Melgiri\\ME559\\Abaqus_CAE\\"
os.chdir(directory)
project_json = json.load(open("C:\\Users\\Nikhil Melgiri\\ME559\\project.json"))


# --[CREATING SKETCH]--
def pin_sketch(model_name, part_name, fraction):
    mdb.Model(name=model_name, modelType=aq.STANDARD_EXPLICIT)
    s1 = mdb.models[model_name].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Line(point1=(-15.0, 0.0), point2=(-20.0, 0.0))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-20.0, 0.0), point2=(-20.0, 20.0))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s1.Line(point1=(-20.0, 20.0), point2=(20.0, 20.0))
    s1.HorizontalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.Line(point1=(20.0, 20.0), point2=(20.0, 0.0))
    s1.VerticalConstraint(entity=g[5], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
    s1.Line(point1=(20.0, 0.0), point2=(15.0, 0.0))
    s1.HorizontalConstraint(entity=g[6], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
    s1.Arc3Points(point1=(-15.0, 0.0), point2=(0.0, 15.0), point3=(-10.0, 11.25))
    s1.Arc3Points(point1=(0.0, 15.0), point2=(15.0, 0.0), point3=(11.25, 10.0))
    s1.EqualLengthConstraint(entity1=g[2], entity2=g[6])
    s1.EqualLengthConstraint(entity1=g[3], entity2=g[5])
    s1.EqualRadiusConstraint(entity1=g[7], entity2=g[8])
    s1.FixedConstraint(entity=g[4])
    s1.undo()
    s1.PerpendicularConstraint(entity1=g[7], entity2=g[2])
    s1.PerpendicularConstraint(entity1=g[8], entity2=g[6])
    s1.ConcentricConstraint(entity1=g[7], entity2=g[8])
    s1.ObliqueDimension(vertex1=v[2], vertex2=v[3], textPoint=(0.294265747070312, 
        26.0482177734375), value=40.0)
    s1.ObliqueDimension(vertex1=v[3], vertex2=v[4], textPoint=(27.954833984375, 
        10.4926605224609), value=20.0)
    s1.RadialDimension(curve=g[8], textPoint=(25.6007461547852, 24.1404571533203), 
        radius=15.0)
    s=mdb.models[model_name].sketches['__profile__']
    s.Parameter(name='Width', expression='40')
    s.Parameter(name='Fraction', expression=str(fraction), previousParameter='Width')
    s.Parameter(name='Length', path='dimensions[0]', expression='Width', 
        previousParameter='Fraction')
    s.Parameter(name='HalfWidth', path='dimensions[1]', expression='Width / 2', 
        previousParameter='Length')
    s.Parameter(name='Radius', path='dimensions[2]', expression='Width * Fraction / 2', 
        previousParameter='HalfWidth')
    p = mdb.models[model_name].Part(name=part_name, dimensionality=TWO_D_PLANAR, 
        type=DEFORMABLE_BODY)
    p = mdb.models[model_name].parts[part_name]
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()
    p = mdb.models[model_name].parts[part_name]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[model_name].sketches['__profile__']


def no_pin_sketch(model_name, part_name, fraction):
    mdb.Model(name=model_name, modelType=aq.STANDARD_EXPLICIT)
    mdb.models[model_name].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    # Create base sketch
    mdb.models[model_name].sketches['__profile__'].Line(point1=(0.0, 15.0), point2=(0.0, 20.0))
    mdb.models[model_name].sketches['__profile__'].VerticalConstraint(addUndoState=False, entity=mdb.models[model_name].sketches['__profile__'].geometry[2])
    mdb.models[model_name].sketches['__profile__'].Line(point1=(0.0, 20.0), point2=(20.0, 20.0))
    mdb.models[model_name].sketches['__profile__'].HorizontalConstraint(addUndoState=False, entity=mdb.models[model_name].sketches['__profile__'].geometry[3])
    mdb.models[model_name].sketches['__profile__'].PerpendicularConstraint(addUndoState=False, entity1=mdb.models[model_name].sketches['__profile__'].geometry[2], entity2=mdb.models[model_name].sketches['__profile__'].geometry[3])
    mdb.models[model_name].sketches['__profile__'].Line(point1=(20.0, 20.0), point2=(20.0, 0.0))
    mdb.models[model_name].sketches['__profile__'].VerticalConstraint(addUndoState=False, entity=mdb.models[model_name].sketches['__profile__'].geometry[4])
    mdb.models[model_name].sketches['__profile__'].PerpendicularConstraint(addUndoState=False, entity1=mdb.models[model_name].sketches['__profile__'].geometry[3], entity2=mdb.models[model_name].sketches['__profile__'].geometry[4])
    mdb.models[model_name].sketches['__profile__'].Line(point1=(20.0, 0.0), point2=(15.0, 0.0))
    mdb.models[model_name].sketches['__profile__'].HorizontalConstraint(addUndoState=False, entity=mdb.models[model_name].sketches['__profile__'].geometry[5])
    mdb.models[model_name].sketches['__profile__'].PerpendicularConstraint(addUndoState=False, entity1=mdb.models[model_name].sketches['__profile__'].geometry[4], entity2=mdb.models[model_name].sketches['__profile__'].geometry[5])
    mdb.models[model_name].sketches['__profile__'].Arc3Points(point1=(0.0, 15.0), point2=(15.0, 0.0), point3=(12.5, 8.75))
    mdb.models[model_name].sketches['__profile__'].EqualLengthConstraint(entity1=mdb.models[model_name].sketches['__profile__'].geometry[3], entity2=mdb.models[model_name].sketches['__profile__'].geometry[4])
    mdb.models[model_name].sketches['__profile__'].EqualLengthConstraint(entity1=mdb.models[model_name].sketches['__profile__'].geometry[2], entity2=mdb.models[model_name].sketches['__profile__'].geometry[5])
    mdb.models[model_name].sketches['__profile__'].ObliqueDimension(textPoint=(16.3181610107422, 24.810962677002), value=20.0, vertex1=mdb.models[model_name].sketches['__profile__'].vertices[1], vertex2=mdb.models[model_name].sketches['__profile__'].vertices[2])
    mdb.models[model_name].sketches['__profile__'].RadialDimension(curve=mdb.models[model_name].sketches['__profile__'].geometry[6], radius=14.3885805415267, textPoint=(25.7375869750977, 21.8998069763184))
    # add parameters
    mdb.models[model_name].sketches['__profile__'].Parameter(expression='20', name='HalfWidth', path='dimensions[0]')
    mdb.models[model_name].sketches['__profile__'].Parameter(expression=str(fraction), name='Fraction', previousParameter='HalfWidth')
    mdb.models[model_name].sketches['__profile__'].Parameter(expression='HalfWidth*Fraction', name='Radius', path='dimensions[1]', previousParameter='Fraction')
    mdb.models[model_name].Part(dimensionality=TWO_D_PLANAR, name=part_name, type=DEFORMABLE_BODY)
    mdb.models[model_name].parts[part_name].BaseShell(sketch=mdb.models[model_name].sketches['__profile__'])
    del mdb.models[model_name].sketches['__profile__']

# --[IMPLEMENTING MATERIAL]--
def implement_material(model_name, part_name, material_name, youngs_modulus, poisson_ratio, density, set_name, section_name):
    mdb.models[model_name].Material(name=material_name)
    mdb.models[model_name].materials[material_name].Density(table=((density, ), ))
    mdb.models[model_name].materials[material_name].Elastic(table=((youngs_modulus, poisson_ratio), ))
    mdb.models[model_name].HomogeneousSolidSection(material=material_name, name=section_name, thickness=None)
    mdb.models[model_name].parts[part_name].Set(faces=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(('[#1 ]', ), ), name=set_name)
    mdb.models[model_name].parts[part_name].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models[model_name].parts[part_name].sets[set_name], sectionName=section_name, thicknessAssignment=FROM_SECTION)
    mdb.models[model_name].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models[model_name].rootAssembly.Instance(dependent=ON, name=part_name, part=mdb.models[model_name].parts[part_name])

# --[IMPLEMENTING STEP]--
def _step(model_name, timestep_name):
    mdb.models[model_name].StaticStep(name=timestep_name, previous='Initial')

# --[IMPLEMENTING BOUNDARY CONDITIONS AND LOADS]--
def pin_implement_bc_loads(model_name, timestep_name, part_name, pressure_magnitude):
    a = mdb.models[model_name].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, predefinedFields=ON, connectors=ON, optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, predefinedFields=OFF, connectors=OFF)
    a = mdb.models[model_name].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models[model_name].parts[part_name]
    a.Instance(name=part_name, part=p, dependent=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, predefinedFields=ON, connectors=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, predefinedFields=OFF, connectors=OFF, adaptiveMeshConstraints=ON)
    mdb.models[model_name].StaticStep(name=timestep_name, previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step=timestep_name)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
    a = mdb.models[model_name].rootAssembly
    s1 = a.instances[part_name].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
    region = regionToolset.Region(side1Edges=side1Edges1)
    mdb.models[model_name].Pressure(name='AppliedLoad', createStepName=timestep_name, region=region, distributionType=UNIFORM, field='', magnitude=-pressure_magnitude, amplitude=UNSET)
    a = mdb.models[model_name].rootAssembly
    e1 = a.instances[part_name].edges
    edges1 = e1.getSequenceFromMask(mask=('[#44 ]', ), )
    region = regionToolset.Region(edges=edges1)
    mdb.models[model_name].YsymmBC(name='BC-1', createStepName=timestep_name, region=region, localCsys=None)
    a = mdb.models[model_name].rootAssembly
    e1 = a.instances[part_name].edges
    edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(edges=edges1)
    mdb.models[model_name].EncastreBC(name='BC-2', createStepName=timestep_name, region=region, localCsys=None)


def no_pin_implement_bc_loads(model_name, timestep_name, part_name, pressure_magnitude):
    mdb.models[model_name].rootAssembly.Surface(name='Surf-1', side1Edges=mdb.models[model_name].rootAssembly.instances[part_name].edges.getSequenceFromMask(('[#4 ]', ), ))
    mdb.models[model_name].Pressure(amplitude=UNSET, createStepName=timestep_name, distributionType=UNIFORM, field='', magnitude=-pressure_magnitude, name='AppliedLoad', region=mdb.models[model_name].rootAssembly.surfaces['Surf-1'])
    mdb.models[model_name].rootAssembly.Set(edges=mdb.models[model_name].rootAssembly.instances[part_name].edges.getSequenceFromMask(('[#10 ]', ), ), name='Set-1')
    mdb.models[model_name].XsymmBC(createStepName=timestep_name, localCsys=None, name='BC-1', region=mdb.models[model_name].rootAssembly.sets['Set-1'])
    mdb.models[model_name].rootAssembly.Set(edges=mdb.models[model_name].rootAssembly.instances[part_name].edges.getSequenceFromMask(('[#2 ]', ), ), name='Set-2')
    mdb.models[model_name].YsymmBC(createStepName=timestep_name, localCsys=None, name='BC-2', region=mdb.models[model_name].rootAssembly.sets['Set-2'])

# --[IMPLEMENTING MESH]--
def pin_mesh(model_name, part_name, top_seeds, bottom_seeds, left_right_seeds, circular_seeds):
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, engineeringFeatures=OFF, mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(meshTechnique=ON)
    
    # top seeds
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    pickedCenterEdges = e.getSequenceFromMask(mask=('[#10 ]', ), )
    p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=pickedCenterEdges, ratio=5.0, number=top_seeds, constraint=FINER)
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#10 ]', ), )
    p.Set(edges=edges, name='TopSeeds')
    
    # left and right seeds
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#28 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=left_right_seeds, constraint=FINER)
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#28 ]', ), )
    p.Set(edges=edges, name='LeftRightSeeds')
    
    # circular seeds
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    pickedEdges1 = e.getSequenceFromMask(mask=('[#2 ]', ), )
    pickedEdges2 = e.getSequenceFromMask(mask=('[#1 ]', ), )
    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
        end2Edges=pickedEdges2, ratio=5.0, number=circular_seeds, constraint=FINER)
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#3 ]', ), )
    p.Set(edges=edges, name='CircularSeeds')
    
    # bottom seeds
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#44 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=bottom_seeds, constraint=FINER)
    p = mdb.models[model_name].parts[part_name]
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#44 ]', ), )
    p.Set(edges=edges, name='BottomSeeds')
    
    # meshing
    mdb.models[model_name].parts[part_name].setMeshControls(elemShape=TRI, regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(('[#1 ]', ), ))
    mdb.models[model_name].parts[part_name].generateMesh()
    mdb.models[model_name].rootAssembly.regenerate()


def no_pin_mesh(top_seeds, bottom_seeds, right_seeds, left_seeds, circular_seeds, model_name, part_name):
    # top seeding
    mdb.models[model_name].parts[part_name].seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end1Edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#8 ]', ), ), number=top_seeds, ratio=5.0)
    mdb.models[model_name].parts[part_name].Set(edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#8 ]', ), ), name='TopSeeding')

    # left seeding
    mdb.models[model_name].parts[part_name].seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end1Edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#10 ]', ), ), number=left_seeds, ratio=5.0)
    mdb.models[model_name].parts[part_name].Set(edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#10 ]', ), ), name='LeftEdgeSeeding')

    # bottom seeding
    mdb.models[model_name].parts[part_name].seedEdgeByNumber(constraint=FINER, edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#2 ]', ), ), number=bottom_seeds)
    mdb.models[model_name].parts[part_name].Set(edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#2 ]', ), ), name='BottomSeeding')

    # right seeding
    mdb.models[model_name].parts[part_name].seedEdgeByNumber(constraint=FINER, edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#4 ]', ), ), number=right_seeds)
    mdb.models[model_name].parts[part_name].Set(edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#4 ]', ), ), name='RightSideSeeding')

    # circular seeding
    mdb.models[model_name].parts[part_name].seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end1Edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#1 ]', ), ), number=circular_seeds, ratio=5.0)
    mdb.models[model_name].parts[part_name].Set(edges=mdb.models[model_name].parts[part_name].edges.getSequenceFromMask(('[#1 ]', ), ), name='CircularSeeding')

    # meshing
    mdb.models[model_name].parts[part_name].setMeshControls(elemShape=TRI, regions=mdb.models[model_name].parts[part_name].faces.getSequenceFromMask(('[#1 ]', ), ))
    mdb.models[model_name].parts[part_name].generateMesh()
    mdb.models[model_name].rootAssembly.regenerate()

# --[IMPLEMENTING AND RUNNING JOB]--
def run_job(job_name, model_name):
    mdb.Job(atTime=None, contactPrint=OFF, description=None, echoPrint=OFF, explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model=model_name, modelPrint=OFF, name=job_name, nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()

# --[POST-PROCESSING]--
def post_process(directory, job_name, model_name, string_identifier, fraction, mesh_level, category, index):
    model = mdb.models[model_name].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=model)
    session.mdbData.summary()
    print(join(directory, '{}.odb'.format(job_name)))
    odb = session.openOdb(name=join(directory, '{}.odb'.format(job_name)))
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S11'), )
    session.printToFile(fileName=join(join(directory, "images"), '{}.png'.format(string_identifier)), format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
    # save s11 data
    data_name = "S11" if index == 1 else "S11_{}".format(index - 1)
    session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
        variable=(('S', INTEGRATION_POINT, ((COMPONENT, "S11"), )), ), 
        operator=MAXIMUM_ENVELOPE, elementSets=(" ALL ELEMENTS", ))
    data = session.xyDataObjects["MAXIMUM_S:{}".format(data_name)]
    part = mdb.models[model_name].parts[part_name]
    max_stress = data.data[-1][1]
    number_of_elements = len(part.elements)
    
    with open(join(join(directory, "data"), '{}.csv'.format(string_identifier)), 'w') as f:
        f.write('category, fraction, mesh level, max stress, number of elements')
        f.write('\n')
        f.write('{},{},{},{},{}'.format(category, fraction, mesh_level, max_stress, number_of_elements))

# --[MAIN]--    
pin = project_json["Pin"]
no_pin = project_json["No Pin"]
pressure = project_json["Pressure"]
fractions = project_json["Fractions"]

model_num = 0
for mesh_data in pin:
    mesh_levels = list(pin["Meshing Levels"].keys())
    for mesh_level in mesh_levels:
        top_seeds = pin["Meshing Levels"][mesh_level]["top_seeds"]
        bottom_seeds = pin["Meshing Levels"][mesh_level]["bottom_seeds"]
        left_right_seeds = pin["Meshing Levels"][mesh_level]["left_right_seeds"]
        for fraction in fractions.keys():
            circular_seeds = fractions[fraction]
            fraction = float(fraction)
            model_num += 1
            string_identifier = "Pin_{}".format(str(model_num))
            model_name = "Model-{}".format(str(model_num))
            job_name = "Pin_Job_{}".format(str(model_num))
            part_name = "Pin_Part_{}".format(str(model_num))
            set_name = "Pin_Set_{}".format(str(model_num))
            timestep_name = "Pin_Timestep_{}".format(str(model_num))
            section_name = "Pin_Section_{}".format(str(model_num))
            pin_sketch(
                model_name=model_name, 
                part_name=part_name, 
                fraction=fraction
            )
            implement_material(
                model_name=model_name,
                part_name=part_name,
                material_name="Isotropic Aluminum",
                youngs_modulus=70e9,
                density=2700,
                set_name=set_name,
                poisson_ratio=0.33,
                section_name=section_name
            )
            _step(model_name=model_name, timestep_name=timestep_name)
            pin_implement_bc_loads(
                model_name=model_name,
                timestep_name=timestep_name,
                part_name=part_name,
                pressure_magnitude=pressure,
            )
            pin_mesh(
                model_name=model_name,
                part_name=part_name,
                top_seeds=top_seeds,
                bottom_seeds=bottom_seeds,
                left_right_seeds=left_right_seeds,
                circular_seeds=circular_seeds
            )
            run_job(job_name=job_name, model_name=model_name)
            post_process(
                directory=directory,
                job_name=job_name,
                model_name=model_name,
                string_identifier=string_identifier,
                fraction=fraction,
                mesh_level=mesh_level,
                category="Pin",
                index=model_num
            )

for mesh_data in no_pin:
    mesh_levels = list(no_pin["Meshing Levels"].keys())
    for mesh_level in mesh_levels:
        top_seeds = no_pin["Meshing Levels"][mesh_level]["top_seeds"]
        bottom_seeds = no_pin["Meshing Levels"][mesh_level]["bottom_seeds"]
        left_seeds = no_pin["Meshing Levels"][mesh_level]["left_seeds"]
        right_seeds = no_pin["Meshing Levels"][mesh_level]["right_seeds"]
        for fraction in fractions:
            circular_seeds = fractions[fraction]
            fraction = float(fraction)
            model_num += 1
            string_identifier = "No_Pin_{}".format(str(model_num))
            model_name = "Model-{}".format(str(model_num))
            job_name = "No_Pin_Job_{}".format(str(model_num))
            part_name = "No_Pin_Part_{}".format(str(model_num))
            set_name = "No_Pin_Set_{}".format(str(model_num))
            timestep_name = "No_Pin_Timestep_{}".format(str(model_num))
            section_name = "No_Pin_Section_{}".format(str(model_num))
            no_pin_sketch(
                model_name=model_name, 
                part_name=part_name, 
                fraction=fraction
            )
            implement_material(
                model_name=model_name,
                part_name=part_name,
                material_name="Isotropic Aluminum",
                youngs_modulus=70e9,
                density=2700,
                set_name=set_name,
                poisson_ratio=0.33,
                section_name=section_name
            )
            _step(model_name=model_name, timestep_name=timestep_name)
            no_pin_implement_bc_loads(
                model_name=model_name,
                timestep_name=timestep_name,
                part_name=part_name,
                pressure_magnitude=pressure,
            )
            no_pin_mesh(
                model_name=model_name,
                part_name=part_name,
                top_seeds=top_seeds,
                bottom_seeds=bottom_seeds,
                right_seeds=right_seeds,
                left_seeds=left_seeds,
                circular_seeds=circular_seeds
            )
            run_job(job_name=job_name, model_name=model_name)
            post_process(
                directory=directory,
                job_name=job_name,
                model_name=model_name,
                string_identifier=string_identifier,
                fraction=fraction,
                mesh_level=mesh_level,
                category="No Pin",
                index=model_num
            )


