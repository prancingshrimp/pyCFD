try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

######################################################################################################################################################################

# for item in xrange(0,1):
for item in xrange(0,129):
    # item = 128
    test_0_vtk = LegacyVTKReader( FileNames=['./test/VTK/test_' + str(item) + '.vtk'] )

    RenderView1 = GetRenderView()
    DataRepresentation2 = Show()
    #DataRepresentation2.EdgeColor = [0.0, 0.0, 0.51]
    #DataRepresentation2.SelectionPointFieldDataArrayName = 'p'
    #DataRepresentation2.SelectionCellFieldDataArrayName = 'cellID'
    #DataRepresentation2.ScalarOpacityUnitDistance = 0.012485316802385881
    #DataRepresentation2.ScaleFactor = 0.034999999403953555

    # a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, RGBPoints=[-6e-6, 0.0, 0.0, 1.0, 3e-5, 1.0, 0.0, 0.0] )

    a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, VectorMode='Magnitude', RGBPoints=[-3e-06, 0.0, 0.0, 1.0, 3e-05, 1.0, 0.0, 0.0] )
    # a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=1, VectorMode='Component', RGBPoints=[-0.0015, 0.0, 0.0, 1.0, 0.0015, 1.0, 0.0, 0.0] )
    # a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=0, VectorMode='Component', RGBPoints=[0.0, 0.0, 0.0, 1.0, 0.051, 1.0, 0.0, 0.0] )

    #RenderView1.CameraPosition = [0.17499999701976776, 0.0, 0.6957348438659209]
    #RenderView1.CameraFocalPoint = [0.17499999701976776, 0.0, 0.0]
    #RenderView1.CameraClippingRange = [0.6290774967616607, 0.7813208648441763]
    #RenderView1.CameraParallelScale = 0.580069427933929
    RenderView1.ViewSize = [1024,768]
    RenderView1.Background = [1.0, 1.0, 1.0]
    # RenderView1.Background = [0.32941176470588235, 0.34901960784313724, 0.42745098039215684]
    RenderView1.CenterAxesVisibility = 0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.CameraParallelProjection = 0
    #RenderView1.CameraFocalPoint = [0.22845810843765293, -0.13974137896955943, 0.0]
    RenderView1.CameraPosition = [0.35/2, 0.0, 0.55]
    RenderView1.CameraFocalPoint = [0.35/2, 0.0, 0.0]


    #DataRepresentation2.ScalarOpacityFunction = []
    DataRepresentation2.ColorArrayName = ('CELL_DATA', 'p')
    DataRepresentation2.LookupTable = a3_U_PVLookupTable
    DataRepresentation2.ColorAttributeType = 'CELL_DATA'

    WriteImage('./pics_for_presentation/p0/p0_' + str(item) + '.png', view=RenderView1)

######################################################################################################################################################################

for item in xrange(0,129):
    # item = 128
    test_0_vtk = LegacyVTKReader( FileNames=['./test/VTK/test_' + str(item) + '.vtk'] )

    RenderView1 = GetRenderView()
    DataRepresentation2 = Show()
    #DataRepresentation2.EdgeColor = [0.0, 0.0, 0.51]
    #DataRepresentation2.SelectionPointFieldDataArrayName = 'p'
    #DataRepresentation2.SelectionCellFieldDataArrayName = 'cellID'
    #DataRepresentation2.ScalarOpacityUnitDistance = 0.012485316802385881
    #DataRepresentation2.ScaleFactor = 0.034999999403953555

    # a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, RGBPoints=[-6e-6, 0.0, 0.0, 1.0, 3e-5, 1.0, 0.0, 0.0] )

    # a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, VectorMode='Magnitude', RGBPoints=[-3e-06, 0.0, 0.0, 1.0, 3e-05, 1.0, 0.0, 0.0] )
    # a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=1, VectorMode='Component', RGBPoints=[-0.0015, 0.0, 0.0, 1.0, 0.0015, 1.0, 0.0, 0.0] )
    a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=0, VectorMode='Component', RGBPoints=[0.0, 0.0, 0.0, 1.0, 0.051, 1.0, 0.0, 0.0] )

    #RenderView1.CameraPosition = [0.17499999701976776, 0.0, 0.6957348438659209]
    #RenderView1.CameraFocalPoint = [0.17499999701976776, 0.0, 0.0]
    #RenderView1.CameraClippingRange = [0.6290774967616607, 0.7813208648441763]
    #RenderView1.CameraParallelScale = 0.580069427933929
    RenderView1.ViewSize = [1024,768]
    RenderView1.Background = [1.0, 1.0, 1.0]
    # RenderView1.Background = [0.32941176470588235, 0.34901960784313724, 0.42745098039215684]
    RenderView1.CenterAxesVisibility = 0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.CameraParallelProjection = 0
    #RenderView1.CameraFocalPoint = [0.22845810843765293, -0.13974137896955943, 0.0]
    RenderView1.CameraPosition = [0.35/2, 0.0, 0.55]
    RenderView1.CameraFocalPoint = [0.35/2, 0.0, 0.0]


    #DataRepresentation2.ScalarOpacityFunction = []
    DataRepresentation2.ColorArrayName = ('CELL_DATA', 'U')
    DataRepresentation2.LookupTable = a3_U_PVLookupTable
    DataRepresentation2.ColorAttributeType = 'CELL_DATA'

    WriteImage('./pics_for_presentation/u0/u0_' + str(item) + '.png', view=RenderView1)

######################################################################################################################################################################

for item in xrange(0,129):
    # item = 128
    test_0_vtk = LegacyVTKReader( FileNames=['./test/VTK/test_' + str(item) + '.vtk'] )

    RenderView1 = GetRenderView()
    DataRepresentation2 = Show()
    #DataRepresentation2.EdgeColor = [0.0, 0.0, 0.51]
    #DataRepresentation2.SelectionPointFieldDataArrayName = 'p'
    #DataRepresentation2.SelectionCellFieldDataArrayName = 'cellID'
    #DataRepresentation2.ScalarOpacityUnitDistance = 0.012485316802385881
    #DataRepresentation2.ScaleFactor = 0.034999999403953555

    # a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, RGBPoints=[-6e-6, 0.0, 0.0, 1.0, 3e-5, 1.0, 0.0, 0.0] )

    # a3_U_PVLookupTable = GetLookupTableForArray( "p", 1, VectorMode='Magnitude', RGBPoints=[-3e-06, 0.0, 0.0, 1.0, 3e-05, 1.0, 0.0, 0.0] )
    a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=1, VectorMode='Component', RGBPoints=[-0.0015, 0.0, 0.0, 1.0, 0.0015, 1.0, 0.0, 0.0] )
    # a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, VectorComponent=0, VectorMode='Component', RGBPoints=[0.0, 0.0, 0.0, 1.0, 0.051, 1.0, 0.0, 0.0] )

    #RenderView1.CameraPosition = [0.17499999701976776, 0.0, 0.6957348438659209]
    #RenderView1.CameraFocalPoint = [0.17499999701976776, 0.0, 0.0]
    #RenderView1.CameraClippingRange = [0.6290774967616607, 0.7813208648441763]
    #RenderView1.CameraParallelScale = 0.580069427933929
    RenderView1.ViewSize = [1024,768]
    RenderView1.Background = [1.0, 1.0, 1.0]
    # RenderView1.Background = [0.32941176470588235, 0.34901960784313724, 0.42745098039215684]
    RenderView1.CenterAxesVisibility = 0
    RenderView1.OrientationAxesVisibility = 0
    RenderView1.CameraParallelProjection = 0
    #RenderView1.CameraFocalPoint = [0.22845810843765293, -0.13974137896955943, 0.0]
    RenderView1.CameraPosition = [0.35/2, 0.0, 0.55]
    RenderView1.CameraFocalPoint = [0.35/2, 0.0, 0.0]


    #DataRepresentation2.ScalarOpacityFunction = []
    DataRepresentation2.ColorArrayName = ('CELL_DATA', 'U')
    DataRepresentation2.LookupTable = a3_U_PVLookupTable
    DataRepresentation2.ColorAttributeType = 'CELL_DATA'

    WriteImage('./pics_for_presentation/v0/v0_' + str(item) + '.png', view=RenderView1)
