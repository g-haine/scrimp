# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/published/MCSS_2023/PV_trace_experiment_012.py
- authors:          Ghislain Haine
- date:             22 sep. 2023
- brief:            trace ParaView post-processings of shallow_water experiments
"""

def trace_experiment(experiment=0, formulation="grad"):
    """
    Post-processing of the three rectangular tank shallow water equations experiments using ParaView
    
    Example of use: to run post-processing of experiment 1 with `grad-grad` formulation, run:
        `python PV_trace_experiment.py 1 grad`
    
    Args:
        - experiment (int): the experiment to reproduce:
            * 0: inviscid case with homogeneous boundary control
            * 1: viscous case with homogeneous boundary control
            * 2: emptying tank by the right end, i.e. normal control prescribed
            * 3: rotating tank, i.e. tangent control prescribed
        - formulation (string): the formulation, i.e., where the first integration by parts is done
            * grad: the first integration by parts is done on the mass preservation equation
            * div: the first integration by parts is done on the linear momentum equation
    """

    if experiment not in [0,1,2,3]:
        raise ValueError(f'Unknown experiment: {experiment}')
    if formulation not in ["grad","div"]:
        raise ValueError(f'Unknown formulation: {formulation}')
    
    # set paths
    import os
    path = os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)),os.path.join("outputs",formulation)),f'experiment_{experiment}',"pv")
    if not os.path.isdir(path):
        raise Exception(f'Experiment {experiment} has not been solved yet.')
    path_h = os.path.join(os.path.join(path,'h'), 'h.pvd')
    path_e_p = os.path.join(os.path.join(path,'e_p'), 'e_p.pvd')

    # trace generated using paraview version 5.10.0-RC1
    #import paraview
    #paraview.compatibility.major = 5
    #paraview.compatibility.minor = 10

    #### import the simple module from the paraview
    import paraview.simple as pvs
    pvs.Connect()
    #### disable automatic camera reset on 'Show'
    pvs._DisableFirstRenderCameraReset()

    # create a new 'PVD Reader'
    hpvd = pvs.PVDReader(registrationName='h.pvd', FileName=path_h)
    hpvd.PointArrays = ['h']

    # get animation scene
    animationScene1 = pvs.GetAnimationScene()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = pvs.GetActiveViewOrCreate('RenderView')

    # reset view to fit data
    renderView1.ResetCamera(True)
    
    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'

    # show data in view
    hpvdDisplay = pvs.Show(hpvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    hpvdDisplay.Representation = 'Surface'
    hpvdDisplay.ColorArrayName = [None, '']
    hpvdDisplay.SelectTCoordArray = 'None'
    hpvdDisplay.SelectNormalArray = 'None'
    hpvdDisplay.SelectTangentArray = 'None'
    hpvdDisplay.OSPRayScaleArray = 'h'
    hpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    hpvdDisplay.SelectOrientationVectors = 'None'
    hpvdDisplay.ScaleFactor = 0.2
    hpvdDisplay.SelectScaleArray = 'None'
    hpvdDisplay.GlyphType = 'Arrow'
    hpvdDisplay.GlyphTableIndexArray = 'None'
    hpvdDisplay.GaussianRadius = 0.01
    hpvdDisplay.SetScaleArray = ['POINTS', 'h']
    hpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    hpvdDisplay.OpacityArray = ['POINTS', 'h']
    hpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    hpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    hpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    hpvdDisplay.ScalarOpacityUnitDistance = 0.2
    hpvdDisplay.OpacityArrayName = ['POINTS', 'h']

    # create a new 'Warp By Scalar'
    warpByScalar1 = pvs.WarpByScalar(registrationName='WarpByScalar1', Input=hpvd)
    warpByScalar1.Scalars = ['POINTS', 'h']

    # show data in view
    warpByScalar1Display = pvs.Show(warpByScalar1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    warpByScalar1Display.Representation = 'Surface'
    warpByScalar1Display.ColorArrayName = [None, '']
    warpByScalar1Display.SelectTCoordArray = 'None'
    warpByScalar1Display.SelectNormalArray = 'None'
    warpByScalar1Display.SelectTangentArray = 'None'
    warpByScalar1Display.OSPRayScaleArray = 'h'
    warpByScalar1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    warpByScalar1Display.SelectOrientationVectors = 'None'
    warpByScalar1Display.ScaleFactor = 0.5
    warpByScalar1Display.SelectScaleArray = 'None'
    warpByScalar1Display.GlyphType = 'Arrow'
    warpByScalar1Display.GlyphTableIndexArray = 'None'
    warpByScalar1Display.GaussianRadius = 0.025
    warpByScalar1Display.SetScaleArray = ['POINTS', 'h']
    warpByScalar1Display.ScaleTransferFunction = 'PiecewiseFunction'
    warpByScalar1Display.OpacityArray = ['POINTS', 'h']
    warpByScalar1Display.OpacityTransferFunction = 'PiecewiseFunction'
    warpByScalar1Display.DataAxesGrid = 'GridAxesRepresentation'
    warpByScalar1Display.PolarAxes = 'PolarAxesRepresentation'
    warpByScalar1Display.ScalarOpacityUnitDistance = 0.5
    warpByScalar1Display.OpacityArrayName = ['POINTS', 'h']

    # hide data in view
    pvs.Hide(hpvd, renderView1)

    # create a new 'Transform'
    transform1 = pvs.Transform(registrationName='Transform1', Input=warpByScalar1)
    transform1.Transform = 'Transform'

    # Properties modified on transform1.Transform
    if experiment in [0,1,2]:
        transform1.Transform.Translate = [0.0, -0.6, -2.]
    else:
        transform1.Transform.Translate = [-2.2, 0.0, -2.]
    transform1.Transform.Scale = [1.0, 1.0, 0.05]

    # show data in view
    transform1Display = pvs.Show(transform1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    transform1Display.Representation = 'Surface'
    transform1Display.ColorArrayName = [None, '']
    transform1Display.SelectTCoordArray = 'None'
    transform1Display.SelectNormalArray = 'None'
    transform1Display.SelectTangentArray = 'None'
    transform1Display.OSPRayScaleArray = 'h'
    transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    transform1Display.SelectOrientationVectors = 'None'
    transform1Display.ScaleFactor = 0.2
    transform1Display.SelectScaleArray = 'None'
    transform1Display.GlyphType = 'Arrow'
    transform1Display.GlyphTableIndexArray = 'None'
    transform1Display.GaussianRadius = 0.01
    transform1Display.SetScaleArray = ['POINTS', 'h']
    transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
    transform1Display.OpacityArray = ['POINTS', 'h']
    transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
    transform1Display.DataAxesGrid = 'GridAxesRepresentation'
    transform1Display.PolarAxes = 'PolarAxesRepresentation'
    transform1Display.ScalarOpacityUnitDistance = 0.2
    transform1Display.OpacityArrayName = ['POINTS', 'h']

    # hide data in view
    pvs.Hide(warpByScalar1, renderView1)

    # set scalar coloring
    pvs.ColorBy(transform1Display, ('POINTS', 'h'))

    # rescale color and/or opacity maps used to include current data range
    transform1Display.UpdatePipeline()
    transform1Display.RescaleTransferFunctionToDataRangeOverTime()

    # show color bar/color legend
    transform1Display.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'h'
    hLUT = pvs.GetColorTransferFunction('h')

    # get color legend/bar for hLUT in view renderView1
    hLUTColorBar = pvs.GetScalarBar(hLUT, renderView1)

    # create a new 'PVD Reader'
    e_ppvd = pvs.PVDReader(registrationName='e_p.pvd', FileName=path_e_p)
    e_ppvd.PointArrays = ['e_p']

    # show data in view
    e_ppvdDisplay = pvs.Show(e_ppvd, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    e_ppvdDisplay.Representation = 'Wireframe'
    e_ppvdDisplay.ColorArrayName = [None, '']
    e_ppvdDisplay.SelectTCoordArray = 'None'
    e_ppvdDisplay.SelectNormalArray = 'None'
    e_ppvdDisplay.SelectTangentArray = 'None'
    e_ppvdDisplay.OSPRayScaleArray = 'e_p'
    e_ppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    e_ppvdDisplay.SelectOrientationVectors = 'None'
    e_ppvdDisplay.ScaleFactor = 0.2
    e_ppvdDisplay.SelectScaleArray = 'None'
    e_ppvdDisplay.GlyphType = 'Arrow'
    e_ppvdDisplay.GlyphTableIndexArray = 'None'
    e_ppvdDisplay.GaussianRadius = 0.01
    e_ppvdDisplay.SetScaleArray = ['POINTS', 'e_p']
    e_ppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    e_ppvdDisplay.OpacityArray = ['POINTS', 'e_p']
    e_ppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    e_ppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    e_ppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    e_ppvdDisplay.ScalarOpacityUnitDistance = 0.2
    e_ppvdDisplay.OpacityArrayName = ['POINTS', 'e_p']

    # create a new 'Glyph'
    glyph1 = pvs.Glyph(registrationName='Glyph1', Input=e_ppvd, GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'e_p']
    glyph1.ScaleArray = ['POINTS', 'e_p']
    if experiment in [0,1,2]:
        glyph1.ScaleFactor = 4.0
    else:
        glyph1.ScaleFactor = 1.0
    glyph1.GlyphTransform = 'Transform2'
    glyph1.GlyphMode = 'All Points'
    
    # show data in view
    glyph1Display = pvs.Show(glyph1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.SelectTCoordArray = 'None'
    glyph1Display.SelectNormalArray = 'None'
    glyph1Display.SelectTangentArray = 'None'
    glyph1Display.OSPRayScaleArray = 'e_p'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.2
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 0.01
    glyph1Display.SetScaleArray = ['POINTS', 'e_p']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'e_p']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'

    # set scalar coloring
    pvs.ColorBy(glyph1Display, ('POINTS', 'e_p', 'Magnitude'))

    # rescale color and/or opacity maps used to include current data range
    glyph1Display.UpdatePipeline()
    glyph1Display.RescaleTransferFunctionToDataRangeOverTime()

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'glyph1'
    e_pLUT = pvs.GetColorTransferFunction('e_p')

    # get color legend/bar for e_pLUT in view renderView1
    e_pLUTColorBar = pvs.GetScalarBar(e_pLUT, renderView1)
    
    if experiment in [0,1,2]:
        # AxesGrid property provides access to the AxesGrid object.
        axesGrid1 = renderView1.AxesGrid

        # To toggle visibility of the axes grid,
        axesGrid1.Visibility = 1
        axesGrid1.XTitle = "Length (m)"
        axesGrid1.XAxisLabels = [0., 0.7, 1.]
        axesGrid1.YTitle = "Width (m)"
        axesGrid1.YAxisLabels = []
        
        # Properties modified on renderView1.AxesGrid
        renderView1.AxesGrid.UseCustomBounds = 1
        renderView1.AxesGrid.CustomBounds = [0.0, 2.0, 0.0, 0.5, 0.0, 0.0]

    # get layout
    layout1 = pvs.GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1612, 811)

    # update the view to ensure updated data information
    renderView1.Update()

    # reset view to fit data
    renderView1.ResetCamera(True)
    if experiment in [0,1,2]:
        renderView1.CameraParallelScale = 0.65
    else:
        renderView1.CameraParallelScale = 1.35
    
    # change scalar bar placement
    # hLUTColorBar.WindowLocation = 'Any Location'
    # hLUTColorBar.Orientation = 'Horizontal'
    # hLUTColorBar.Position = [0.1, 0.03]
    # hLUTColorBar.ScalarBarLength = 0.33
    
    # change scalar bar placement
    # e_pLUTColorBar.WindowLocation = 'Any Location'
    # e_pLUTColorBar.Orientation = 'Horizontal'
    # e_pLUTColorBar.Position = [0.55, 0.03]
    # e_pLUTColorBar.ScalarBarLength = 0.33

    # get the time-keeper
    timeKeeper1 = pvs.GetTimeKeeper()

    for k in range(0, len(timeKeeper1.TimestepValues), 25):
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timeKeeper1.TimestepValues[k]

        # update the view to ensure updated data information
        renderView1.Update()

        # save screenshot
        pvs.SaveScreenshot(os.path.join(path, 'time_'+str(timeKeeper1.TimestepValues[k])+'.png'), renderView1, ImageResolution=[1612, 811], 
            # PNG options
            CompressionLevel='0')
                
    if experiment in [0,1,2]:
        # hide axes grid
        axesGrid1.Visibility = 0

    # create a new 'Compute Derivatives'
    computeDerivatives1 = pvs.ComputeDerivatives(registrationName='ComputeDerivatives1', Input=e_ppvd)
    computeDerivatives1.Scalars = [None, '']
    computeDerivatives1.Vectors = ['POINTS', 'e_p']
    computeDerivatives1.OutputVectorType = 'Vorticity'
    computeDerivatives1.OutputTensorType = 'Nothing'

    # show data in view
    computeDerivatives1Display = pvs.Show(computeDerivatives1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    computeDerivatives1Display.Representation = 'Surface'
    computeDerivatives1Display.ColorArrayName = [None, '']
    computeDerivatives1Display.SelectTCoordArray = 'None'
    computeDerivatives1Display.SelectNormalArray = 'None'
    computeDerivatives1Display.SelectTangentArray = 'None'
    computeDerivatives1Display.OSPRayScaleArray = 'e_p'
    computeDerivatives1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    computeDerivatives1Display.SelectOrientationVectors = 'Vorticity'
    computeDerivatives1Display.ScaleFactor = 0.2
    computeDerivatives1Display.SelectScaleArray = 'None'
    computeDerivatives1Display.GlyphType = 'Arrow'
    computeDerivatives1Display.GlyphTableIndexArray = 'None'
    computeDerivatives1Display.GaussianRadius = 0.01
    computeDerivatives1Display.SetScaleArray = ['POINTS', 'e_p']
    computeDerivatives1Display.ScaleTransferFunction = 'PiecewiseFunction'
    computeDerivatives1Display.OpacityArray = ['POINTS', 'e_p']
    computeDerivatives1Display.OpacityTransferFunction = 'PiecewiseFunction'
    computeDerivatives1Display.DataAxesGrid = 'GridAxesRepresentation'
    computeDerivatives1Display.PolarAxes = 'PolarAxesRepresentation'
    computeDerivatives1Display.ScalarOpacityUnitDistance = 0.20854610816029576
    computeDerivatives1Display.OpacityArrayName = ['POINTS', 'e_p']
    computeDerivatives1Display.SelectInputVectors = ['POINTS', 'e_p']
    computeDerivatives1Display.WriteLog = ''

    # hide data in view
    pvs.Hide(e_ppvd, renderView1)

    # hide data in view
    pvs.Hide(glyph1, renderView1)

    # hide data in view
    pvs.Hide(transform1, renderView1)

    # set scalar coloring
    pvs.ColorBy(computeDerivatives1Display, ('CELLS', 'Vorticity', 'Z'))

    # rescale color and/or opacity maps used to include current data range
    computeDerivatives1Display.UpdatePipeline()
    computeDerivatives1Display.RescaleTransferFunctionToDataRangeOverTime()

    # show color bar/color legend
    computeDerivatives1Display.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'Vorticity'
    VorticityLUT = pvs.GetColorTransferFunction('Vorticity')

    # get color legend/bar for VorticityLUT in view renderView1
    VorticityLUTColorBar = pvs.GetScalarBar(VorticityLUT, renderView1)
    
    # change scalar bar placement
    VorticityLUTColorBar.WindowLocation = 'Any Location'
    VorticityLUTColorBar.Orientation = 'Horizontal'
    VorticityLUTColorBar.Position = [0.33, 0.03]
    VorticityLUTColorBar.ScalarBarLength = 0.33

    # update the view to ensure updated data information
    renderView1.Update()

    # reset view to fit data
    renderView1.ResetCamera(True)
    if experiment==3:
        renderView1.CameraParallelScale = 1.5

    for k in range(0, len(timeKeeper1.TimestepValues), 25):
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timeKeeper1.TimestepValues[k]

        # update the view to ensure updated data information
        renderView1.Update()
        
        # save screenshot
        pvs.SaveScreenshot(os.path.join(path, 'vorticity_time_'+str(timeKeeper1.TimestepValues[k])+'.png'), renderView1, ImageResolution=[1612, 811], 
            # PNG options
            CompressionLevel='0')
                
    # 3D Movie
    # Properties modified on transform1.Transform
    transform1.Transform.Translate = [0.0, 0.0, -2.]
    transform1.Transform.Scale = [1.0, 1.0, 0.05]

    # show data in view
    pvs.Show(e_ppvd, renderView1)
    pvs.Show(glyph1, renderView1)
    pvs.Show(transform1, renderView1)
    
    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)
    transform1Display.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    pvs.Hide(computeDerivatives1, renderView1)    

    #change interaction mode for render view
    renderView1.InteractionMode = '3D'

    # get layout
    layout1 = pvs.GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1544, 848)

    # current camera placement for renderView1
    if experiment in [0,1,2]:
        renderView1.CameraPosition = [2.1263496101156196, -4.4520034029147135, 1.3514417760354152]
        renderView1.CameraFocalPoint = [1.0000030994406826, 0.25000328128589905, 0.24985101123456846]
        renderView1.CameraViewUp = [-0.10467515299388723, 0.20301325361172984, 0.9735649599300942]
        renderView1.CameraViewAngle = 16.335814722911497
        renderView1.CameraParallelScale = 1.060716148402028
    else:
        renderView1.CameraPosition = [1.3695991534047924, -5.716939728390033, 1.6499235130602161]
        renderView1.CameraFocalPoint = [8.982419967651367e-05, 0.0001659393310546875, 0.3105143569409847]
        renderView1.CameraViewUp = [-0.10467515299388723, 0.20301325361172984, 0.9735649599300942]
        renderView1.CameraViewAngle = 17.09823263460748
        renderView1.CameraParallelScale = 1.560549425066
        

    # update the view to ensure updated data information
    renderView1.Update()

    animationScene1.GoToFirst()

    #animationScene1.Play()

    pvs.SaveAnimation(os.path.join(path,'movie.ogv'), pvs.GetActiveView(), ImageResolution=[1544, 848], FrameRate=20)

    renderView1.ResetActiveCameraToNegativeZ()
    animationScene1.GoToFirst()
    
    # close everything
    pvs.RemoveViewsAndLayouts()
    pvs.ResetSession()
    pvs.Disconnect()

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        experiment = int(sys.argv[1])
        if len(sys.argv)>2:
            formulation = sys.argv[2]
    else:
        experiment = 0
        formulation = "grad"
    trace_experiment(experiment, formulation)

