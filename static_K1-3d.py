import os
import numpy
import sys
import platform
from operator import *
from abaqus import *
from abaqusConstants import *
from part import *
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
from odbAccess import *
import math

# input params


L   = float(8)
H   = float(14)
t   = float(1)
rho = float(1000)
E   = float(200)
nu  = float(0.3)
F   = float(2)
cl  = float(2)


model_name = 'toy_model'

min_r = 1.e-3
r = 0.01


Ncontour = 5
Nhoop = 4
MeshCircle = 0.01
MinEdgeMesh = H*0.1
MaxEdgeMesh = H*0.1
epsilon = 0.02
def calculate_arc_midpoint(center, start_point, end_point):
    """
    Calculate the midpoint of a circular arc.

    Parameters:
        center (tuple): Coordinates of the arc's center (cx, cy, cz).
        start_point (tuple): Coordinates of the arc's start point (sx, sy, sz).
        end_point (tuple): Coordinates of the arc's end point (ex, ey, ez).

    Returns:
        tuple: Coordinates of the midpoint of the arc (mx, my, mz).
    """
    # Extract coordinates
    cx, cy, cz = center
    sx, sy, sz = start_point
    ex, ey, ez = end_point

    # Calculate the radius of the arc
    radius = math.sqrt((sx - cx)**2 + (sy - cy)**2 + (sz - cz)**2)

    # Calculate the angles of the start and end points relative to the center
    start_angle = math.atan2(sy - cy, sx - cx)
    end_angle = math.atan2(ey - cy, ex - cx)

    # Calculate the midpoint angle (average of start and end angles)
    mid_angle = (start_angle + end_angle) / 2

    # Calculate the midpoint coordinates
    midpoint = (
        cx + radius * math.cos(mid_angle),
        cy + radius * math.sin(mid_angle),
        cz  # Assuming the arc lies in a plane, z-coordinate remains constant
    )

    return midpoint
mid1=calculate_arc_midpoint((L-cl,H/2,0),(L-cl-cl/2,H/2,0),calculate_arc_midpoint((L-cl,H/2,0),(L-cl,H/2-cl/2,0),(L-cl-cl/2,H/2,0)))
# Ref: http://bertoldi.seas.harvard.edu/files/bertoldi/files/learnabaqusscriptinonehour.pdf
print(mid1)
model = mdb.Model(name=model_name)
part = model.Part(name='Part 1', dimensionality=THREE_D, 
                  type=DEFORMABLE_BODY)

box_sketch = model.ConstrainedSketch(name='cantilever', sheetSize=1)
box_sketch.Line(point1=(0,0),point2=(L,0))
box_sketch.Line(point1=(L,0),point2=(L,H))
box_sketch.Line(point1=(L,H),point2=(0,H))
box_sketch.Line(point1=(0,H),point2=(0,0))
part.BaseSolidExtrude(sketch=box_sketch,depth=1)

# Sketch Partitions
my_part = mdb.models['toy_model'].parts['Part 1']
f= my_part.faces
pickedFaces = f.findAt((2,0,0),)
edges=my_part.edges
edge=edges.findAt((5.,0.,0.),)
cell = my_part.cells[0]
##Open sheet
partition_sketch = model.ConstrainedSketch(name='Partition Sketch', sheetSize=1)
##Sketching crack
partition_sketch.Line(point1=(L,H/2), point2=(L-cl,H/2))
##Sketching Circles
partition_sketch.CircleByCenterPerimeter(center=(L-cl,H/2), point1=(L-cl+r,H/2))
partition_sketch.CircleByCenterPerimeter(center=(L-cl,H/2), point1=(L-cl+(cl/2),H/2))
##Sketching Square lines
partition_sketch.Line(point1=(L,(H/2)-cl),point2=(L,+(H/2)+cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)+cl),point2=(L,(H/2)+cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)-cl),point2=(L-2*cl,(H/2)+cl))
####Sketching 90 degree lines
partition_sketch.Line(point1=(L-2*cl,(H/2)-cl),point2=(L,(H/2)-cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)),point2=(L-cl,(H/2)))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-cl,(H/2)+cl))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-cl,(H/2)-cl))
#####Sketching 45 degree lines
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L,(H/2)+cl))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L,(H/2)-cl))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-2*cl,(H/2)+cl))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-2*cl,(H/2)-cl))
my_part.PartitionFaceBySketch(faces=pickedFaces, sketch=partition_sketch)
e1=my_part.edges
##3d partition
edges=(e1.findAt((L-(cl/4),H/2.,0.),),e1.findAt((L-cl+r/2,H/2.,0.),),e1.findAt((L-(3*cl/4),H/2.,0.),))
smallcircle=(e1.findAt((L-cl,H/2+r/2,0.),),e1.findAt((L-cl-r/2,H/2+r/2.,0.),),e1.findAt((L-cl+r/2,H/2+r/2.,0.),),e1.findAt((L-cl+r/2,H/2-r/2.,0.),),e1.findAt((L-cl-r/2,H/2-r/2.,0.),),e1.findAt((L-cl,H/2-r/2.,0.),))

my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=cell,edges=edges,sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=cell,edges=(e1.findAt((L-cl,H/2+r/2,0.),),e1.findAt((L-cl-r/2,H/2+r/2.,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=cell,edges=(e1.findAt((L-cl+r/2,H/2+r/2.,0.)),e1.findAt((L-cl+r/2,H/2-r/2.,0.),),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=cell,edges=(e1.findAt((L-cl-r/2,H/2-r/2.,0.),),e1.findAt((L-cl,H/2-r/2.,0.),)),sense=FORWARD)


my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=cell,edges=e1.findAt((L-cl,H/2+r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl,H/2+r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl-r,H/2+r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl+r,H/2+r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl+r,H/2-r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl-r,H/2-r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=e1.findAt((L-cl,H/2-r,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=(e1.findAt((L-cl,H/2+r,0.),),e1.findAt((L-cl-r,H/2+r,0.),)),sense=FORWARD)

my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl,H/2+cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl-cl/2,H/2+cl/2,0.),),e1.findAt((L-cl-cl/2,H/2+cl,0.),),e1.findAt((L-cl,H/2+cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl,H/2-cl/2,0.),),e1.findAt((L-cl-cl/2,H/2-cl,0.),),e1.findAt((L-2*cl,H/2-cl/2,0.),),e1.findAt((L-2*cl,H/2+cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl,H/2-cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl,H/2-cl,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl/2,H/2+cl,0.),),e1.findAt((L-cl/2,H/2+cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl/2,H/2,0.),),e1.findAt((L-r,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl/2,H/2-cl/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl/2,H/2,0.),),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-3*cl/2,H/2-cl/2,0.),),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl,H/2+cl/2,0.),),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl-cl/2,H/2,0.),),e1.findAt((L-cl-cl/2+r,H/2,0.),),e1.findAt((L-cl-r/2,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=e1.findAt((L-cl-cl/2,H/2,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=e1.findAt((L-cl,H/2-cl/2,0.),),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl+r,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl+r,H/2,0.),),e1.findAt((L-cl,H/2-r,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[7],edges=(e1.findAt((L-cl,H/2-r,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[6],edges=(e1.findAt((L-cl-r,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[6],edges=(e1.findAt((L-cl-cl/2,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells,edges=(e1.findAt((L-cl-r,H/2,0.),)),sense=FORWARD)
my_part.PartitionCellByExtrudeEdge(line=e1.findAt((L,0.,500.E-03),),cells=my_part.cells[1],edges=(e1.findAt((L-cl-r,H/2,0.),)),sense=FORWARD)


##material and section
material = model.Material(name='steel')
material.Density(table=((rho, ), ))
material.Elastic(table=((E, nu), ))
model.HomogeneousSolidSection(material='steel', name='steel_section')

## Section Assignments
region = my_part.Set(cells=my_part.cells, name='All_Cells')
my_part.SectionAssignment(region=region, sectionName='steel_section',
 offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)



## Instance
assembly = mdb.models['toy_model'].rootAssembly
instance = assembly.Instance(name='Instance 1', part=part, dependent=OFF)


# Convert the datum plane to cell partition
datum_plane = assembly.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=L/2)
datum_id = datum_plane.id
partition_face = assembly.datums[datum_id]
assembly.PartitionCellByDatumPlane(datumPlane=partition_face, cells=instance.cells)

##Step
step_name = 'Step-1'
mdb.models['toy_model'].StaticStep(name=step_name, previous='Initial', timePeriod=1.0)
# Output
model.FieldOutputRequest(
name='F-Output-1',
createStepName='Step-1',
variables=('S', 'U', 'RF', 'ENER', 'ELEN'))
##History and field output
model.fieldOutputRequests['F-Output-1'].setValues(variables=('S','U','RF', 
                                                                        'ENER','ELEN'))
model.HistoryOutputRequest(contourType=J_INTEGRAL, contourIntegral='Crack-tip', 
                                createStepName='Step-1', name='J-Integral', 
                                numberOfContours=Ncontour, rebar=EXCLUDE, 
                                sectionPoints=DEFAULT, 
                                stressInitializationStep='Initial');
model.HistoryOutputRequest(contourType=K_FACTORS, contourIntegral='Crack-tip', 
                                createStepName='Step-1', name='K-value', 
                                numberOfContours=Ncontour, rebar=EXCLUDE, 
                                sectionPoints=DEFAULT, 
                                stressInitializationStep='Initial');
model.HistoryOutputRequest(contourType=T_STRESS, contourIntegral='Crack-tip', 
                                createStepName='Step-1', name='T-Stress', 
                                numberOfContours=Ncontour, rebar=EXCLUDE, 
                                sectionPoints=DEFAULT, 
                                stressInitializationStep='Initial')


##crack
face_indexes=[22,42,79]
f1 = assembly.instances['Instance 1'].faces
picked_faces = instance.faces[22:23] + instance.faces[42:43] + instance.faces[79:80]
assembly.Set(faces=picked_faces, name='CRACK_SEAM')
pickedRegions = assembly.sets['CRACK_SEAM']
mdb.models['toy_model'].rootAssembly.engineeringFeatures.assignSeam(
    regions=pickedRegions)
crack_vertices = instance.vertices.findAt(((L-cl,H/2,0),))
crack_edges = instance.edges.findAt(((6.,7.,500.E-03),))
##Defining Crack with contour integral
assembly.engineeringFeatures.ContourIntegral(
    collapsedElementAtTip=SINGLE_NODE,
    extensionDirectionMethod=Q_VECTORS,
    symmetric=OFF,
    midNodePosition=0.25,
    crackFront=Region(edges=crack_edges),
    crackTip=Region(edges=crack_edges),
    name=('Crack-tip'),
    qVectors=(((L,H/2,0.0), (L-cl,H/2,0.0)),)
    )

##Couple
##reference point
ref_point_coords = (L/2,H,500.E-03)
ref_point =assembly.ReferencePoint(point=ref_point_coords)
ref_point_id = ref_point.id
# Create assembly set for the reference point
assembly.Set(referencePoints=(assembly.referencePoints[ref_point_id],), name='RP_Set')
##edge for couple
facecouple1 = instance.faces.findAt(((2*L/3,H,0.7),))  
facecouple2 = instance.faces.findAt(((L/6,H,0.7),)) 
facecouple=(facecouple1,facecouple2) 
assembly.Set(faces=facecouple, name='face_Set')


# Create the coupling constraint
# Access the model
mdb.models['toy_model'].Coupling(
    name='Coupling_Constraint',
    controlPoint=assembly.sets['RP_Set'],
    surface= assembly.sets['face_Set'],
    influenceRadius=WHOLE_SURFACE,
    couplingType=KINEMATIC,
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON  # Adjust degrees of freedom as needed
)


# fixed BC
BC1 = instance.faces.findAt(((2*L/3,0,0.7),))
BC2 = instance.faces.findAt(((L/6,0,0.7),))
facesbc=(BC1,BC2)
region = assembly.Set(faces=facesbc, name='fixedBC_set')
fixedBC = model.DisplacementBC(name='fixed_BC',createStepName='Initial',
region=region, u1=0.0, u2=0.0, u3=UNSET)



##cells
c1 = assembly.instances['Instance 1'].cells
cell_indices = [2, 6,7,8,13,24,25,26]
pickedcell =   tuple(instance.cells[i] for i in cell_indices)

# Meshing Parameters
assembly.seedPartInstance(regions=(instance, ), deviationFactor=0.1, 
                          minSizeFactor=0.1, size=0.2)
assembly.setMeshControls(regions=pickedcell, technique=SWEEP,elemShape=WEDGE)
cell_indices2 = [0,1,3,4,5,9,10,11,12,14,15,16,17,18,19,20,21,22,23]
pickedcell2 =   tuple(instance.cells[i] for i in cell_indices2)

assembly.setMeshControls(regions=pickedcell2 , elemShape=HEX, 
                         technique=STRUCTURED)


# Assign Element Type
assembly.setElementType(elemTypes=(ElemType(elemCode=C3D8,elemLibrary=STANDARD),
                                   ElemType(elemCode=C3D8, 
                                            elemLibrary=STANDARD)), 
                        regions=(instance.cells, ))


# Select edges for local seeding
# Replace the coordinates with the actual location of the edges
edges_to_seed = instance.edges.findAt(
    ((0, H, 0.5),),  # First edge coordinates
    ((L, H, 0.5),),
    ((0, 0, 0.5),),
    ((L, 0, 0.5),),  # Second edge coordinates (if applicable)
)

# Assign local seed by number to the selected edges
number_of_elements = 6  # Desired number of elements along the edge
assembly.seedEdgeByNumber(edges=edges_to_seed, number=number_of_elements, constraint=FINER)


# Generate the mesh
assembly.generateMesh(regions=instance.cells,seedConstraintOverride=ON,
                      meshTechniqueOverride=OFF)


#Snpashot
session.viewports['Viewport: 1'].setValues(displayedObject=assembly)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON,
                                                     optimizationTasks=OFF, 
                                                     geometricRestrictions=OFF, 
                                                     stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(meshTechnique=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=30.,
                                                farPlane=31,
                                                width=L*2,
                                                height=L*2,
                                                viewOffsetX=0,
                                                viewOffsetY=0)
session.printToFile(fileName='mesh_01.png', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))

load_values = [7000, 200, 300]
with open('results.dat', 'w') as f:
    for i, load_value in enumerate(load_values):

        # Apply the concentrated force
        mdb.models['toy_model'].ConcentratedForce(
            name='Load-{}'.format(i + 1), createStepName=step_name, region=assembly.sets['RP_Set'],
            cf1=load_value, distributionType=UNIFORM, field='', localCsys=None
        )
        mdb.models['toy_model'].Moment(name='LoadMoment-{}'.format(i + 1),
        createStepName=step_name, region=assembly.sets['RP_Set'],cm3=141055)
        
        # Create a job for the current load
        job_name = 'CrackJob-{}'.format(i + 1)
        mdb.Job(name=job_name, model=model_name, type=ANALYSIS, explicitPrecision=SINGLE,
                nodalOutputPrecision=SINGLE, description='CrackJob-{}'.format(i + 1))
        
        # Submit the job
        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()

        
        # Extracting results
        session.openOdb(name='CrackJob-{}'.format(i + 1)+'.odb',  readOnly=True)
        output = session.odbs['CrackJob-{}'.format(i + 1)+'.odb']

        # Access the step in the ODB
        step = output.steps[step_name]

        for historyRegionName in step.historyRegions.keys():
            print(historyRegionName)


        # To find keys, 1) abaqus cae nogui 2) >>> print(output.steps.keys()) 
        # 3) >>> print(output.steps['Step-1'].historyRegion.keys()) ...
        tmp = output.steps['Step-1'].historyRegions['ElementSet  ALL ELEMENTS'].historyOutputs

        K1_value = []
        K2_value = []
        J_value = []
        T_value = []


        for item in tmp.keys():
            if item[0:2] == 'K1':
                K1_value.append(tmp[item].data[0][1])   
            if item[0:2] == 'K2':
                K2_value.append(tmp[item].data[0][1])    
            if item[0:2] == 'J ':
                J_value.append(tmp[item].data[0][1])   
            if item[0:2] == 'T-':
                T_value.append(tmp[item].data[0][1])

        load=load_values[i]
            # Write headers
        f.write('# K1 value CrackJob-{}\n'.format(load))
        for item in K1_value:
            f.write('%.4e ' % item)
        f.write('\n')

        f.write('# K2 value CrackJob-{}\n'.format(load))
        for item in K2_value:
            f.write('%.4e ' % item )
        f.write('\n')

        f.write('# J value CrackJob-{}\n'.format(load))
        for item in J_value:
            f.write('%.4e ' % item)
        f.write('\n')

        f.write('# T value CrackJob-{}\n'.format(load))
        for item in T_value:
            f.write('%.4e ' % item)
        f.write('\n')