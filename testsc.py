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
print(edge)
partition_sketch = model.ConstrainedSketch(name='Partition Sketch', sheetSize=1)
partition_sketch.Line(point1=(L,H/2), point2=(L-cl,H/2))
partition_sketch.CircleByCenterPerimeter(center=(L-cl,H/2), point1=(L-cl+r,H/2))
partition_sketch.CircleByCenterPerimeter(center=(L-cl,H/2), point1=(L-cl+(cl/2),H/2))
partition_sketch.Line(point1=(L,(H/2)-cl),point2=(L,+(H/2)+cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)+cl),point2=(L,(H/2)+cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)-cl),point2=(L-2*cl,(H/2)+cl))
####
partition_sketch.Line(point1=(L-2*cl,(H/2)-cl),point2=(L,(H/2)-cl))
partition_sketch.Line(point1=(L-2*cl,(H/2)),point2=(L-cl,(H/2)))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-cl,(H/2)+cl))
partition_sketch.Line(point1=(L-cl,(H/2)),point2=(L-cl,(H/2)-cl))
#####
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
my_part.SectionAssignment(region=region, sectionName='steel_section', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)



## Instance
assembly = mdb.models['toy_model'].rootAssembly
instance = assembly.Instance(name='Instance 1', part=part, dependent=OFF)


# Convert the datum plane to cell partition
datum_plane = assembly.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=L/2)
datum_id = datum_plane.id
partition_face = assembly.datums[datum_id]
assembly.PartitionCellByDatumPlane(datumPlane=partition_face, cells=instance.cells)


# Step
model.StaticStep(name='Step-1', previous='Initial')

##crack
f1 = assembly.instances['Instance 1'].faces
pickedface=instance.faces.getSequenceFromMask(mask=('[#0:2 #8000 ]',))+instance.faces.getSequenceFromMask(mask=('[#0 #200 ]',),)+instance.faces.getSequenceFromMask(mask=('[#0 #200 ]',),)+instance.faces.getSequenceFromMask(mask=('[#0 #400 ]',),)

assembly.Set(faces=pickedface, name='CRACK_SEAM')
pickedRegions = assembly.sets['CRACK_SEAM']
mdb.models['toy_model'].rootAssembly.engineeringFeatures.assignSeam(regions=pickedRegions)
crack_vertices = instance.vertices.findAt(((L-cl,H/2,0),))
assembly.engineeringFeatures.ContourIntegral(
    collapsedElementAtTip=SINGLE_NODE,
    extensionDirectionMethod=Q_VECTORS,
    symmetric=OFF,
    midNodePosition=0.25,
    crackFront=Region(vertices=crack_vertices),
    crackTip=Region(vertices=crack_vertices),
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
edgecouple = instance.edges.findAt(((L/2,H,0.7),))  # Replace with coordinates on the edge
assembly.Set(edges=(edgecouple,), name='Edge_Set')


# Create the coupling constraint
# Access the model
mdb.models['toy_model'].Coupling(
    name='Coupling_Constraint',
    controlPoint=assembly.sets['RP_Set'],
    surface= assembly.sets['Edge_Set'],
    influenceRadius=WHOLE_SURFACE,
    couplingType=KINEMATIC,
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON  # Adjust degrees of freedom as needed
)
##Loading
mdb.models['toy_model'].ConcentratedForce(name='Load-1', createStepName='Step-1', 
                        region=assembly.sets['RP_Set'], cf2=-F, distributionType=UNIFORM, 
                        field='', localCsys=None)

# fixed BC
edgeBC1 = instance.edges.findAt(((L,0.,500.E-03),))
edgeBC2 = instance.edges.findAt(((0,0.,500.E-03),))
edgesbc=(edgeBC1,edgeBC2)
region = assembly.Set(edges=edgesbc, name='fixedBC_set')
fixedBC = model.DisplacementBC(name='fixed_BC',createStepName='Initial',region=region, u1=0.0, u2=0.0, u3=UNSET)

# Output
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
                           stressInitializationStep='Initial');