

%[baronDuFourierH, baronDuFourierV, xcorFig] = brownianPointingSim(verid_delta_scalar,rot_delta_scalar, fft_horiz_graph, fft_vert_graph, xcorr_graph, nRows, nColumns, plot_num)

master_bDuFH = figure;
master_bDuFV = figure;
master_xcor = figure;
master_xCorrPos = figure;
nRows = 2;
nColumns = 3;

brownianPointingSim(0.7,0,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 1); %verid_under
brownianPointingSim(0,0.7,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 2); %rot_under
brownianPointingSim(0.35,0.35,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 3); %both_under
brownianPointingSim(1.1,0,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 4); %verid_over
brownianPointingSim(0,1.1,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 5); %rot_over
brownianPointingSim(0.55,0.55,master_bDuFH, master_bDuFV, master_xcor,master_xCorrPos, nRows, nColumns, 6); %both_over
