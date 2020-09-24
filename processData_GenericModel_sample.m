function processData_GenericModel_sample

% This function processes the Loughborough dataset gait data using
% traditional or generic modelling methods. The processes included here
% are:
%
%   - Importing and converting data
%   - Scaling the model using generic methods
%   - Inverse kinematics
%   - Residual reduction
%   - Inverse muscle tracking (OpenSim Moco)
%
% Note that this script should be run from it's own directory to ensure the
% relative file paths are appropriate.

%% Set-up

%Import opensim libraries
import org.opensim.modeling.*

%Turn annoying warnings off
warning off

%Add supplementary code folder to MATLAB path. The double dot represents
%the previous folder
addpath(genpath('..\Supplementary'));

%Add model geometry directory. 'pwd' identifies current folder
ModelVisualizer.addDirToGeometrySearchPaths([pwd,'\..\..\Generic\Geometry']);

%% Get and convert data
%  TODO: consider packaging some of the heavier code that will be applied
%  to multiple trials into one or more functions...

%Set participant ID
%TODO: edit this to allow use across multiple participants (e.g. could use
%the input function to request participant ID)
prompt = 'What is the participantID?';
participantID = input(prompt,'s')

%Navigate to participants data directory. 'cd' changes current folder
cd(['..\..\PublicDataset\',participantID,'\Gait\']);
participantDir = [pwd,'\'];

%Load in the static trial and convert marker data to .trc

%Load trial
%TODO: check that this file naming convention persists across all
%participants. 
%In this dataset the static trials fro P01, P02, P03 etc. are
%named P1ANT, P2ANT, P3ANT etc. What this short if statement is doing is
%checking whether the participant number is below 10, and if so is removing
%the leading ‘0’ from the filename. This isn’t needed for participants
%above number 10, as there files are named P10ANT, P11ANT etc.
%if str2num(participantID(2:end)) < 10
%   staticFile = fullfile(pwd,[participantID(1:2:end),'ANT.c3d']);
%else
%   staticFile = fullfile(pwd,[participantID,'ANT.c3d']);
%end
% 1 = COP
staticFile = fullfile(pwd,'WBDS01static1.c3d');
staticC3D = osimC3D(staticFile,1);

%Rotate marker data to align with OpenSim coordinate system
%staticC3D.rotateData('x',-90);
%staticC3D.rotateData('z',-90);

%Convert mm to m units
staticC3D.convertMillimeters2Meters();

%There are markers in the c3d file that we don't need + ones that cause
%errors once converted. We can set a list of markers to keep and then
%delete the unwanted columns
staticKeepMarkers = [{'C7'}; {'T10'}; {'L1'}; {'STERNUM'}; {'XIPHOID'}; 
    {'R.ASIS'}; {'R.PSIS'}; {'L.ASIS'}; {'L.PSIS'};{'R.GTR'}; {'RHJC'}; 
    {'RT1'}; {'RT2'}; {'RT3'}; {'RT4'}; {'LT1'}; {'LT2'}; {'LT3'}; {'LT4'};
    {'RKJC'}; {'RS1'}; {'RS2'}; {'RS3'}; {'RS4'};{'RTOE'};
    {'LHJC'}; {'LKJC'}; {'LS1'}; {'LS2'}; {'LS3'}; {'LS4'};{'LTOE'}
    {'R.Knee'}; {'R.Knee.Medial'}; {'R.Ankle'}; {'R.Ankle.Medial'}; {'R.Heel'};
    {'R.MT1'}; {'R.MT5'}; {'L.GTR'};{'L.Knee'}; {'L.Knee.Medial'}; 
    {'L.Ankle'}; {'L.Ankle.Medial'}; {'L.Heel'}; {'L.MT1'}; {'L.MT5'}; {'L.HF'}; {'L.Iliac.Crest'};
    {'L.MT2'}; {'L.TT'}; {'R.HF'}; {'R.Iliac.Crest'}; {'R.MT2'}; {'R.TT'} ];

%Get the marker table
markerTable = staticC3D.getTable_markers();

%Loop through the list of markers and identify those that need to be removed
staticRemoveList = {};
for mm = 0:staticC3D.getTable_markers().getColumnLabels.size()-1
    %Get current marker name
    currMarker = char(staticC3D.getTable_markers().getColumnLabels.get(mm));
    %Check if current marker isn't in the keep list
    if ~any(strcmp(staticKeepMarkers,currMarker))
        %Add to remove list
        staticRemoveList = [staticRemoveList, currMarker];
    end
end
clear mm
    
%Loop through the delete list and remove them from the marker table
for dd = 1:length(staticRemoveList)
    markerTable.removeColumn(staticRemoveList{dd})    
end
clear dd

%Write marker table to file using TRCFileAdapter. The staticFile variable
%is a string pointing to the full path and name of the static trial file.
%This includes it’s file extension – ‘.c3d’. We want to name this new TRC
%file the same name, but obviously change the extension because it is a
%.trc file – so this is why we take off four letters from the end (i.e. to
%remove the .c3d) and then concatenate our new file extension .trc onto it.
TRCFileAdapter.write(markerTable,[staticFile(1:end-4),'.trc']);


%Set participant mass and height 
%TODO: this will also need to be input somehow as it will vary
prompt = 'What is the mass of participant?';
participantmass = input(prompt,'s')
prompt = 'What is the height of participant?';
participantheight = input(prompt,'s')

%Next we'll convert the dynamic trials to .trc and .mot data so that they
%can be processed later.
%TODO: this script currently only processes one designated trial for P01.
%It will need to be edited so that it takes in the trial names required,
%and loops through these appropriately.
gaitFile = fullfile(pwd,'WBDS01walkT01.c3d');

%Load gait trial
gaitC3D = osimC3D(gaitFile,1);

%Rotate marker and force data to align with OpenSim coordinate system
%gaitC3D.rotateData('x',-90);
%gaitC3D.rotateData('z',-90);

%Convert mm to m units
gaitC3D.convertMillimeters2Meters();

%As before there are markers in the C3D file we might not want. We can
%remove them in a similar fashion to only take what we need. The dynamic
%trials don't have the hip and knee joint centres
gaitKeepMarkers = [{'C7'}; {'T10'}; {'L1'}; {'STERNUM'}; {'XIPHOID'}; 
    {'R.ASIS'}; {'R.PSIS'}; {'L.ASIS'}; {'L.PSIS'};{'R.GTR'}; {'RHJC'}; 
    {'RT1'}; {'RT2'}; {'RT3'}; {'RT4'}; {'LT1'}; {'LT2'}; {'LT3'}; {'LT4'};
    {'RKJC'}; {'RS1'}; {'RS2'}; {'RS3'}; {'RS4'};{'RTOE'};
    {'LHJC'}; {'LKJC'}; {'LS1'}; {'LS2'}; {'LS3'}; {'LS4'};{'LTOE'}
    {'R.Knee'}; {'R.Knee.Medial'}; {'R.Ankle'}; {'R.Ankle.Medial'}; {'R.Heel'};
    {'R.MT1'}; {'R.MT5'}; {'L.GTR'};{'L.Knee'}; {'L.Knee.Medial'}; 
    {'L.Ankle'}; {'L.Ankle.Medial'}; {'L.Heel'}; {'L.MT1'}; {'L.MT5'}; {'L.HF'}; {'L.Iliac.Crest'};
    {'L.MT2'}; {'L.TT'}; {'R.HF'}; {'R.Iliac.Crest'}; {'R.MT2'}; {'R.TT'}];

%Get the marker table
markerTable = gaitC3D.getTable_markers();

%Loop through the list of markers and identify those that need to be removed
gaitRemoveList = {};
for mm = 0:gaitC3D.getTable_markers().getColumnLabels.size()-1
    %Get current marker name
    currMarker = char(gaitC3D.getTable_markers().getColumnLabels.get(mm));
    %Check if current marker isn't in the keep list
    if ~any(strcmp(gaitKeepMarkers,currMarker))
        %Add to remove list
        gaitRemoveList = [gaitRemoveList, currMarker];
    end
end
clear mm
    
%Loop through the delete list and remove them from the marker table
for dd = 1:length(gaitRemoveList)
    markerTable.removeColumn(gaitRemoveList{dd})    
end
clear dd

%Write marker table to file using TRCFileAdapter
TRCFileAdapter.write(markerTable,[gaitFile(1:end-4),'.trc']);

%Next we need to write the GRF data to a .mot and .xml file to be used in
%the simulations. This is relatively easy for a trial like this which has
%one foot contact, as it's easy to allocate the force to the appropriate
%body. It becomes a little more difficult with multiple foot contacts or
%forces (e.g. running on a treadmill).

%Get the forces table from the C3D. This has the same f1, p1, m1 column
%headers as before. We need to convert these to individual columns of
%force_x, _y and _z (and the same for the position and moment columns).
forceTable = gaitC3D.getTable_forces();

%We can flatten the Vec3 table to a TimeSeriesTable. This is an OpenSim
%function (flatten) that converts the Vec3 version of a time series table
%to a standard time series table. The Vec3 tables has each column
%containing three components (e.g. good for marker data). To create a .mot
%file for our GRF data we need a normal time series table, so we ‘flatten’
%it.
forceTimeSeries = forceTable.flatten();

%Create a StdVectorString to fix column labels. StdVectorString is an
%OpenSim class that basically contains a list of strings.
forceLabels = [{'ground_force_vx1'}; {'ground_force_vy1'}; {'ground_force_vz1'};
    {'ground_force_px1'}; {'ground_force_py1'}; {'ground_force_pz1'};
    {'ground_torque_x1'}; {'ground_torque_y1'}; {'ground_torque_z1'};
    {'ground_force_vx2'}; {'ground_force_vy2'}; {'ground_force_vz2'};
    {'ground_force_px2'}; {'ground_force_py2'}; {'ground_force_pz2'};
    {'ground_torque_x2'}; {'ground_torque_y2'}; {'ground_torque_z2'}];
forceString = StdVectorString();
for ff = 1:length(forceLabels)
    forceString.add(forceLabels{ff});
end
clear ff

%Set the labels in the TimeSeriesTable. The variable forceTimeSeries is a TimeSeriesTable, which is an OpenSim class that stores data in a tabular format with a time index associated to it. We set the column labels just so they are consistent with what we want.
forceTimeSeries.setColumnLabels(forceString);

%The centre of pressure and moment data has NaN's instead of zeroes, which
%could be problematic down the line. Run a quick fix for this. This is just
%a convenient way to check. In reality the x, y and z component will have
%nan’s in them – but we only need to check one to understand whether the
%nan’s are present. You could equally check in the y and z and get the same
%result – or even check that each one has a nan to be sure – but I think
%just checking the x is fine in this instance. Note that this was something
%unique I came across about this dataset, I haven’t had to do this before
%with my own data.We can also set force, moment and pressure values to zero
%when the vertical force doesn't meet a specified threshold (> 20N here).
for pp = 0:forceTimeSeries.getNumRows()-1
    %Check for nans in point data
    if isnan(forceTimeSeries.getDependentColumn('ground_force_px1').get(pp))
        %Convert to zeroes
        forceTimeSeries.getDependentColumn('ground_force_px1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_py1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_pz1').set(pp,0);
    end
     if isnan(forceTimeSeries.getDependentColumn('ground_force_px2').get(pp))
        %Convert to zeroes
        forceTimeSeries.getDependentColumn('ground_force_px2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_py2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_pz2').set(pp,0);
    end
    %Check for nans in torque data
    if isnan(forceTimeSeries.getDependentColumn('ground_torque_x1').get(pp))
        %Convert to zeroes
        forceTimeSeries.getDependentColumn('ground_torque_x1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_y1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_z1').set(pp,0);
    end
     if isnan(forceTimeSeries.getDependentColumn('ground_torque_x2').get(pp))
        %Convert to zeroes
        forceTimeSeries.getDependentColumn('ground_torque_x2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_y2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_z2').set(pp,0);
    end
    %Check whether data doesn't meet force threshold
    if forceTimeSeries.getDependentColumn('ground_force_vy1').get(pp) < 20
        %Set all data to zeroes
        forceTimeSeries.getDependentColumn('ground_force_px1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_py1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_pz1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_x1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_y1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_z1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vx1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vy1').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vz1').set(pp,0);
    end
     if forceTimeSeries.getDependentColumn('ground_force_vy2').get(pp) < 20
        %Set all data to zeroes
        forceTimeSeries.getDependentColumn('ground_force_px2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_py2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_pz2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_x2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_y2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_torque_z2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vx2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vy2').set(pp,0);
        forceTimeSeries.getDependentColumn('ground_force_vz2').set(pp,0);
    end
end
clear pp

%Filter force data
%TODO: this just uses a 50Hz low pass Butterworth filter to tidy up the
%noisiness, we should do a better job of selecting the filter

%Get force data into Matlab structure to filter
for ff = 0:forceTimeSeries.getColumnLabels().size()-1
    %Get the current column
    d = forceTimeSeries.getDependentColumn(forceTimeSeries.getColumnLabels().get(ff));
    %Loop through and extract data to a matrix
    for gg = 0:d.size()-1
        forceData(gg+1,ff+1) = d.get(gg);
    end
    clear gg
end
clear ff

%Create butterworth filter
fc = 50;
fs = gaitC3D.getRate_force();
filterOrder = 4;
[b,a] = butter(filterOrder,fc/(fs/2));

%Filter force data. The above filtfilt command is a Matlab function that
%filters (i.e. smooths data). This code is applying the low-pass filter we
%created (allocated to the b and a variables) to our force data to smooth
%it. In this instance it doesn’t look at the second, element, but looks
%across the 2nd dimension (i.e. columns) and gets the size (i.e. the number
%of columns). This loop just then goes through column by column and filters
%the data until it reaches the end.
for ff = 1:size(forceData,2)
    forceData(:,ff) = filtfilt(b,a,forceData(:,ff));
end
clear ff

%Put filtered data back in TimeSeriesTable
for ff = 0:forceTimeSeries.getColumnLabels().size()-1
    %Loop through and put relevant data back into column
    for gg = 0:forceTimeSeries.getDependentColumn(forceTimeSeries.getColumnLabels().get(ff)).size()-1
        forceTimeSeries.getDependentColumn(forceTimeSeries.getColumnLabels().get(ff)).set(gg,forceData(gg+1,ff+1));
    end
    clear gg
end
clear ff

%Write the force data to a .mot file using the STOFileAdapter
STOFileAdapter.write(forceTimeSeries,[gaitFile(1:end-4),'_grf.mot']);

%Next we need to create the external loads XML file associted with the MOT data

%Create a blank external loads object
grfLoads = ExternalLoads();

%Set name of the external loads object. 'trialName' will be the filename
[~, trialName] = fileparts(gaitFile);
grfLoads.setName(trialName);

%Create an external force object representing the GRFs
grfForce = ExternalForce();
%Identify the foot the force needs to go to based on the trial name
%if contains(trialName,'R(')
%    forceBody = 'calcn_r';
%elseif contains(trialName,'L(')
%    forceBody = 'calcn_l';
%end

%Set force name
%grfForce.setName(forceBody);

%Set the body to apply the force to
%grfForce.set_applied_to_body(forceBody);

%Set the force, point and torque string identifiers. These are the
%‘identifiers’ that are present in our GRF mot data that identifies what
%each column refers to. If you check the .mot data file you’ll see that
%there are column labels that link to these identifiers with x, y and z
%added to them. What happens in OpenSim when loading external load data is
%it loads the .mot file, but it doesn’t know what each column of data
%refers to – so it looks to the .xml file we create and sees that the
%columns labelled ‘ground_force_vx’, ‘ground_force_vy’ and
%‘ground_force_vz’ contain the force data – and so on for the other data
%types. 
grfForce.set_force_identifier('ground_force_v');
grfForce.set_point_identifier('ground_force_p');
grfForce.set_torque_identifier('ground_torque_');

%Set data source to 'Unassigned' as the data will come from the MOT file.
%The ‘Unassigned’ aspect is just a default component we add. I think
%we do this because we set the data to come from the .mot file in another
%aspect of the external loads .xml file. I assume that if you set the data
%another way you might set something differently here.
grfForce.set_data_source_name('Unassigned');

%Append the force to the external loads object
grfLoads.cloneAndAppend(grfForce);

%Set the datafile for the external loads to the .mot file created
%Note this designation to the file is just the name, and doesn't include a
%file path, but will still work because they are in the same directory as
%the .mot file.
grfLoads.setDataFileName([trialName,'_grf.mot']);

%Print the external loads to an XML file
grfLoads.print([gaitFile(1:end-4),'_grf.xml']);

% Now from the above we have everything converted into the necessary
% formats to use in OpenSim commands!

%% Model scaling

%Navigate to participants processed directory to store outputs. Somewhere
%previous to this we would have allocated a string that represents a folder
%to the ‘participantDir’ variable (it sounds like this is the P01 folder).
%This ‘cd’ command changes the working directory to the participantDir,
%then goes up a directory, and then goes into the ‘ProcessedGeneric’
%folder. 'mkdir' makes a new folder, 'cd' changes the current folder.
%'processedDir = [pwd,'\']' is just allocating the current directory to
%a variable called ‘processedDir.’ ‘pwd’ is a matlab command for returning the current
%directory as a string. I add the extra ‘\’ as I prefer my directory path
%variables to include a trailing slash – I just think this makes it easier
%to use in later commands.
cd([participantDir,'..\ProcessedGeneric']); processedDir = [pwd,'\'];
mkdir('Scaling_Sample'); cd('Scaling_Sample');

%Set up a scale tool
scaleTool = ScaleTool();

%Set subject mass
massKg=str2num(participantmass);
scaleTool.setSubjectMass(massKg);

%Set generic model file (unscaled model). The four dots and slashes is
%going up four directories from the current directory and then into the
%‘GenericMSM’ folder.
scaleTool.getGenericModelMaker().setModelFileName('..\..\..\..\Generic\LaiArnold2017_refined_torso_Loughborough_PublicDataset.osim');

%Set the measurement set. 'ModelScaler' - specifying parameters for
%sclaling the model. 'MeasurementSet' Specifies the measurements by which
%body segments are to be scaled (measurement-based scaling). 'ScaleSet' for
%manual scaling. These .xml files are available and will be in the
%‘GenericMSM’ folder (you should be able to view them using something like
%Notepad++). Generally you can adapt an existing tasks or measurement set
%from a past project to suit your current model and/or marker set – or if
%you’re using the same model and marker set you can re-use it. Creating
%these files is one of the early jobs to do in running a project.
measurementSet = MeasurementSet.safeDownCast(OpenSimObject.makeObjectFromFile('..\..\..\..\Generic\scaleMeasurementSet.xml'));
scaleTool.getModelScaler().setMeasurementSet(measurementSet);

%Set scale tasks. 'MarkerPlacer' Specifies parameters for placing markers
%on the model once a model is scaled. 'IKTaskSet' Task set used to specify weights used
%in the IK computation of the static pose.
taskSet = IKTaskSet('..\..\..\..\Generic\scaleTasks.xml');
for k = 0:taskSet.getSize()-1
    scaleTool.getMarkerPlacer().getIKTaskSet().cloneAndAppend(taskSet.get(k));
end
clear k

%Set marker file. TRC file (.trc) containing the marker positions is used for
%measurement-based scaling. This is usually a static trial, but doesn't
%need to be.  The marker-pair distances are computed for each time step in
%the TRC file and averaged across the time range.
[~,staticName] = fileparts(staticFile);
scaleTool.getMarkerPlacer().setMarkerFileName(['..\..\Gait\',staticName,'.trc']);
scaleTool.getModelScaler().setMarkerFileName(['..\..\Gait\',staticName,'.trc']);

%Set options. 'ModelScaler' - specifying parameters for scaling the
%model.'ScalingOrder' specifies method (measurements or manualScale) and
%order of scaling. Similar to above, the scale order parameter in the scale
%tool needs to be set using an OpenSim ‘ArrayStr’ object. In this case, we
%have an ‘ArrayStr’ object called scaleOrder. We then set the first string
%(i.e. index 0) in this object to be ‘measurements’ which is the term we
%want to use for scaling order (I think another option here is to use
%manual scaling or something). So the zero here just refers to the index of
%the string, not a numerical value of zero.
scaleTool.getModelScaler().setPreserveMassDist(true);
scaleOrder = ArrayStr(); scaleOrder.set(0,'measurements');
scaleTool.getModelScaler().setScalingOrder(scaleOrder);

%Set time range (0.5 - 1.5 seconds. You would think that there would be a
%command that lets you specify the two times, but you need to use the
%OpenSim ‘ArrayDouble’ class to do this. What these few commands do is
%create an ‘ArrayDouble’ object called timeRange. We then set the first
%time value (i.e. index 0) as 0.5; and then set the second time (i.e. index
%1) as 1.5. This therefore sets a time range from 0.5-1.5.
timeRange = ArrayDouble();
timeRange.set(0,0.5); timeRange.set(1,1.5);
scaleTool.getMarkerPlacer().setTimeRange(timeRange);
scaleTool.getModelScaler().setTimeRange(timeRange);

%Set output files. 'output_model_file' generates the scaled-only OpenSim
%model file after scaling and before markers are moved. 'output_scale_file'
%generates a 'ScaleSet' file consisting of the x-y-z scale factors used to
%scale each segment.'GenericScaledModel.osim' (optional)generates the
%scaled-only OpenSim model file (using 'sclaeMeasurementSet'.xml' )after
%scaling and before markers are moved.'GenericScaleSet.xml' generates a
%'ScaleSet' file (using 'sclaeMeasurementSet'.xml' )consisting of the x-y-z
%scale factors used to scale each
%segment.'GenericScaledModel_adjusted.osim' generates OpenSim model
%file(using 'scaleTasks.xml' ) after scaling and marker
%placement.'StaticMotion.mot' (optional) generates motion file (using
%'scaleTasks.xml') after marker relocation.'StaticMarkers.xml' generates
%output marker set file (using 'scaleTasks.xml') containing the new marker
%locations after markers have been placed.

%This can be useful for understanding/looking back at how the different
%segments have been scaled in different axes and relative to one another.
%Remember I mentioned that it’s good if these are realistic and somewhat
%similar across the different segments.
scaleTool.getModelScaler().setOutputModelFileName('GenericScaledModel_sample.osim');

%This can be useful for understanding/looking back at how the different
%segments have been scaled in different axes and relative to one another.
%Remember I mentioned that it’s good if these are realistic and somewhat
%similar across the different segments.
scaleTool.getModelScaler().setOutputScaleFileName('GenericScaleSet_sample.xml');

%The ‘adjusted’ here refers to the adjustment of the markers we allow to
%move on the model (i.e. those set as ‘false’ for the fixed parameter in
%the marker set). We don’t allow the markers that should be on consistent
%spots (i.e. anatomical landmarks) to move. We do allow things like
%tracking markers on segments that will likely be in different positions
%from participant-to-participant to move to the appropriate spots.
scaleTool.getMarkerPlacer().setOutputModelFileName('GenericScaledModel_adjusted_sample.osim');

%this is like an inverse kinematics output for the static trial. In
%reality, it should just be the person standing still.
scaleTool.getMarkerPlacer().setOutputMotionFileName('StaticMotion_sample.mot');

%this isn’t really that much of a useful output, and wouldn’t be used that
%often. In the ‘_adjusted’ .osim model generated above, it includes the
%marker set in the model file – which is simply what this static markers
%file contains.
scaleTool.getMarkerPlacer().setOutputMarkerFileName('StaticMarkers_sample.xml');

%Save and run scale tool. All the above mentioned output files are
%generated using 'GenericScaleSetup.xml'. When you run the setup file
%through the scale tool (i.e. scaleTool.run()), these outputs are
%generated.
scaleTool.print('GenericScaleSetup_sample.xml');
scaleTool.run();

% Scaling and marker placement looks OK - knee markers aren't placed that
% great with the experimental data, so there placement/weighting could be
% reviewed. 

%Load in scaled model
scaledModel = Model('GenericScaledModel_adjusted_sample.osim');

%Scale muscle strength based on linear function presented in Handsfield et
%al. (2014). This uses some convenience functions that are packaged with
%the Rajagopal et al. (2015) gait model. Note that the height of the
%generic model is 1.700 and the height of the experimental participant is
%set earlier. This is using a supplementary
%function‘scaleOptimalForceSubjectSpecific’ to scale the models maximum
%muscle forces based on the linear relationship between height and mass.
%You can see this function in the Code > Supplementary directory. It is
%comparing the height and mass of the generic model to our scaled model,
%and then scaling the muscle forces in the scaled model accordingly. The
%two ‘finalize’ commands are common things you have to do when you change a
%model like we have and before you print it out. The last command gives us
%an extra model file that includes the changes to the maximum muscle
%forces.
heightM=str2num(participantheight);
genModel = Model('..\..\..\..\Generic\LaiArnold2017_refined_torso_Loughborough_PublicDataset.osim');
scaledModelMuscle = scaleOptimalForceSubjectSpecific(genModel,scaledModel,1.700,heightM);
scaledModelMuscle.finalizeFromProperties();
scaledModelMuscle.finalizeConnections();
scaledModelMuscle.print('GenericScaledModel_StrengthScaled_sample.osim');

%% Inverse kinematics

%Navigate up and create a directory to store IK results
cd('..'); mkdir('IK_sample'); cd('IK_sample'); ikDir = [pwd,'\'];

%Initialise IK tool
ikTool = InverseKinematicsTool();

%Set name
ikTool.setName(trialName);

%Set model
ikTool.setModel(scaledModelMuscle);

%Set task set
ikTool.set_IKTaskSet(IKTaskSet('..\..\..\..\GenericMSM\ikTasks.xml'));

%Set marker file
[~,trialName] = fileparts(gaitFile);
ikTool.set_marker_file(['..\..\Gait\',trialName,'.trc']);

%Set times. This is more just a consideration of Meg’s project. We are
%processing this data to look at tibial stress, which is only really
%relevant during the stance phase. For other projects, you might take a
%whole gait cycle, in which you would have to adjust the event times that
%you’re looking for. For example, if you had someone running on a treadmill
%you would more so look for consecutive foot strikes (rather than foot
%strike and toe off like we have here). In summary, this start and end time
%setting will be dependent on what you’re interested in for a specific
%project. Identify time based on the stance phase of the trial, which we'll
%identify from a threshold from the vertical GRF. These forces can be
%grabbed from the current trials .mot file
gaitGRF = TimeSeriesTable(['..\..\Gait\',trialName,'_grf.mot']);
vGRF = gaitGRF.getDependentColumn('ground_force_vy');
%Loop through and extract data to a matrix for ease of use
for gg = 0:vGRF.size()-1
    vGRF_data(gg+1) = vGRF.get(gg);
end
clear gg

%Find the first instance where the vertical GRF is > 20N
startGRFind = find(vGRF_data > 20,1);

%Find where GRF returns below 50N following the start point
endGRFind = find(vGRF_data(startGRFind:end) < 20,1) + startGRFind;

%Identify the times associated with these indices in the GRF data
startTime = gaitGRF.getIndependentColumn().get(startGRFind-1);
endTime = gaitGRF.getIndependentColumn().get(endGRFind-1);

%Set start and end times in IK tool
ikTool.setStartTime(startTime);
ikTool.setEndTime(endTime);

%Set reporting of marker model positions to compare to experimental data if
%desired. Setting this to true will generate the
%‘_model_marker_locations.sto’ file when the IK tools is run. I believe
%there may be another setting in the IK tool like
%‘set_report_marker_errors’. This may be set to true by default so we
%didn’t need to change it. Something you could investigate perhaps.
ikTool.set_report_marker_locations(true);

%Set IK outputs
ikTool.set_output_motion_file([trialName,'_ik_sample.mot']);

%Run IK
ikTool.print(['ikSetup_sample_',trialName,'.xml']);
clc
ikTool.run();

%TODO: add some form of plotting function here to view results...

%% Residual Reduction

%Before running muscle driven simulations, we should run iterations of the
%residual reduction algorithm to make sure the kinematics are consistent
%with the GRF data.

%NOTE: the actuator and tasks files used here are directly taken from the
%gait2392 model in the OpenSim resources - these could be tweaked a little
%for better performance

%NOTE: subtalar angle is currently unlocked, but mtp angle is locked...

%Navigate and create an RRA directory
cd('..'); mkdir('RRA_sample'); cd('RRA_sample'); rraDir = [pwd,'\'];

%Initialise RRA Tool
rraTool = RRATool();

%Set general settings outside of loop

%Set tool to replace model force set
rraTool.setReplaceForceSet(true);

%Set actuators file. I think this stems from the fact that, if you want,
%you can set multiple force set files. An ArrayStr object allows you to set
%multiple file paths as different indices of the array. So theoretically if
%you had your actuators split across multiple files, you could do this.
%Seemingly this wouldn’t be possible with the tasks file.The first command
%here creates an empty ArrayStr – so it has no ‘strings’ in it. Because
%there’s nothing in it, we set the first string at index 0.
forceSetFiles = ArrayStr();
forceSetFiles.set(0,'..\..\..\..\..\GenericMSM\rraActuators.xml');
rraTool.setForceSetFiles(forceSetFiles);

%Set tracking tasks file
rraTool.setTaskSetFileName('..\..\..\..\..\GenericMSM\rraTasks.xml');

%Set a low pass filter frequency on the kinematics data (12Hz is OK for
%running I think...)
%TODO: consider best filtering practices for this data
rraTool.setLowpassCutoffFrequency(12);

%Set to adjust the COM to reduce residuals
rraTool.setAdjustCOMToReduceResiduals(true);

%Set the torso body COM to be the one that gets adjusted
rraTool.setAdjustedCOMBody('torso');

%Set external loads file
rraTool.setExternalLoadsFileName(['..\..\..\Gait\',trialName,'_grf.xml']);

%Loop through and perform rra three times
%Each time we'll adjust the mass specified by the RRA tool, and then re-use
%this adjusted model in the next iteration. On the first iteration, we'll
%use the IK motion data - but on subsequent iterations we'll use the
%adjusted rra kinematics. Similarly, on the first iteration we'll set the
%model to our scaled model - but then use the adjusted version on
%subsequent iterations.
for rr = 1:3
    
    %Set tool name based on iteration
    rraTool.setName(['rra_sample',num2str(rr)]);
    
    %Go to an appropriate results directory
    mkdir(['rra_sample',num2str(rr)]); cd(['rra_sample',num2str(rr)]);

    %Set desired kinematics file
    if rr == 1
        
        %Use IK data
        rraTool.setDesiredKinematicsFileName(['..\..\IK_sample\',trialName,'_ik_sample.mot']);
        
        %Set initial and final time
        %You need to use the IK first and last times here as I don't think the tool
        %likes if the IK data doesn't have the times you put in
        rraTool.setStartTime(Storage(['..\..\IK_sample\',trialName,'_ik_sample.mot']).getFirstTime());
        rraTool.setFinalTime(Storage(['..\..\IK_sample\',trialName,'_ik_sample.mot']).getLastTime());
        
        %Set to original scaled model
        rraTool.setModelFilename('..\..\Scaling\GenericScaledModel_StrengthScaled_sample.osim');
    
    else
        
        %Use RRA result from last iteration
        %Using the raw file was causing RRA to crash - so we can write some
        %code to load it in and print as a new sto file. A reason for this
        %might have been how large the RRA file was, so we resample here.
        %Get the kinematics as a sotra
        rraKinematics = Storage(['..\rra_sample',num2str(rr-1),'\rra_sample',num2str(rr-1),'_Kinematics_q.sto']);
        
        %Resample the kinematics. So the 250 is kind of arbitrarily decided
        %– it will downsample the kinematic time stamp (RRA uses really
        %small time steps) to a point that still retains smooth kinematics.
        %I think the main reason for this is for the later Moco procedure –
        %the resampling also evenly spaces out the time-steps. In Moco, if
        %there are duplicate (i.e. shared time steps) you get an error.
        %This can sometimes happen with the small time steps RRA produces,
        %so here we can kind of fix this problem before it eventuates.
        rraKinematics.resampleLinear(1/250);
        
        %Ensure the file is set to read degrees, as this is what the output
        %is. Again, I think this is just a check to make sure the data goes
        %to Moco correctly. I think I found that sometimes when you load
        %the RRA kinematics in it assumes that it is in radians, when we
        %had it in degrees. This just makes sure it sees it as degrees
        %which is correct.
        rraKinematics.setInDegrees(true);
        
        %Print to a new file in the current directory
        rraKinematics.print('previousRRA_kinematics_sample.sto');
        
        %Set this in the rra tool
        rraTool.setDesiredKinematicsFileName('previousRRA_kinematics_sample.sto');
        
        %Reset the times to the new file
        rraTool.setStartTime(rraKinematics.getFirstTime());
        rraTool.setFinalTime(rraKinematics.getLastTime());
        
        %Set to adjusted model from previous run
        rraTool.setModelFilename(['..\rra_sample',num2str(rr-1),'\rraAdjustedModel_sample_',num2str(rr-1),'.osim']);
    
    end
    
    %Set output model file
    rraTool.setOutputModelFileName(['rraAdjustedModel_sample_',num2str(rr),'.osim']);
    
    %Print out the tool
    rraTool.print(['setupRRA_sample',num2str(rr),'.xml']);
    
    %Run rra tool
    %Sometimes the typical .run() command doesn't work well. I don't know
    %why... You can use this opensim-cmd command
% % %     rraTool.run();
    clc
    system(['opensim-cmd run-tool setupRRA_sample',num2str(rr),'.xml'])
    
    %TODO: create some code to examine specifics & outputs of each RRA run
    
    %RRA automatically adjusts the torso bodies COM, but it doesn't apply
    %the recommended mass adjustments. These need to be extracted from the
    %tool log and then applied to the adjusted model. This approach again
    %uses a convenience function provided with the Rajagopal (2015) model
    %tools.
    
    %Get the outlog and read the mass and OCM changes. This is just a check
    %of the data in that it displays the changes being made to the model.
    %The first highlighted line reads in the out.log file that is created
    %when RRA is run. The other two highlighted lines are just using a bit
    %of code to piece out the two lines of text that we’re interested –
    %these are consistent in each RRA file in how the out.log displays the
    %change in body mass and the change in the centre of mass position for
    %the adjusted body (the torso in our case).
    outlog = fileread('out.log');
    massChange = str2double(outlog(regexpi(outlog, 'total mass change: ', 'end'):regexpi(outlog, 'total mass change: .?[0-9]+[.][0-9]+', 'end')));
    disp(['RRA_sample',num2str(rr),' dMass = ', num2str(massChange)])
    dCOM = outlog(regexpi(outlog, 'Mass Center \(COM\) adjustment:', 'end'):regexpi(outlog, 'Mass Center \(COM\) adjustment: .+]', 'end'));
    disp(['RRA_sample',num2str(rr),' dCOM = ', dCOM])
    
    %Apply the mass changes to the model using the convenience function
    osimModel_rraMassChanges = Model(rraTool.getModelFilename());
    osimModel_rraMassChanges = setMassOfBodiesUsingRRAMassChange(osimModel_rraMassChanges, massChange);
    osimModel_rraMassChanges.finalizeFromProperties();
    osimModel_rraMassChanges.finalizeConnections();
    osimModel_rraMassChanges.print(rraTool.getOutputModelFileName());
    
    %Navigate back up to general RRA directory
    cd('..');    
    
end
clear rr

% By the end of the 3 RRA iterations we should have relatively low
% residuals, and oberserve only a small mass change required. I did notice
% with this example that the third iteration seemed to introduce some larger
% residuals in places, but improvement in others still. The mass change was
% somewhat larger for the 3rd vs. 2nd iteration, and the torso COM continued
% to progressively move further back (even outside of the torso body). 
%
% TODO: A more ideal approach is to iteratively check your RRA results
% against some sort of threshold that you're happy with and cease the RRA
% iterations at that point (e.g. peak, average residuals, mass change etc.)
