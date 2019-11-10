function mstrials = jnm_elim_ms( session, evts )

% If adding trials that were not complete (i.e. did not achieve juice),
% evts are the list of other event codes to look for. They will be looked
% for in the sequence you enter them. Default is to rfori looking for each
% successfully completed presentation in reverse order. 

if ~exist( 'evts', 'var' )
    evts = 32 : -2 : 24;
end

plotms = 0;
rawpath = jnm_rawpath();

bhv = concatBHV( [ rawpath session '.bhv'] );

folder = '/Users/jakew/Dropbox/data/eyemove/';
samplerate = 1000;

mstrials = nan( 3, length( bhv.AnalogData ) );

for trial = 1 : length( bhv.AnalogData )
    
    disp( ['Begin Detection: Trial no. ' num2str( trial ) ] )
    
    isabort = 0;
    
    samples = cat(2, (1:size(bhv.AnalogData{1,trial}.EyeSignal, 1)).', ...
        bhv.AnalogData{1,trial}.EyeSignal( :, 1:2 ), ...
        nan(size(bhv.AnalogData{1,trial}.EyeSignal, 1),1), ...
        nan(size(bhv.AnalogData{1,trial}.EyeSignal, 1),1));
    
    if size( bhv.AnalogData{1,trial}.EyeSignal, 2 ) > 2
        
        samples( :, 4:5 ) = bhv.AnalogData{1,trial}.EyeSignal( :, 4:5 );
        
    end
    
    if any( bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == 96 ) )
        
        samples = samples( bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == 23 ) - 100 : ...
            bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == 96 ), : );
        
    else
        
        for npres = [ evts, NaN ]
            
            if any( bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == npres ) )
                
                samples = samples( bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == 23 ) - 100 : ...
                    bhv.CodeTimes{ 1, trial }( bhv.CodeNumbers{ 1, trial } == npres ), : );
                
                break
                
            end
            
            if isnan( npres )
                
                isabort = 1;
                
            end
            
        end 
    end
    
    if isabort
        
        disp( ['DID NOT MAKE FIXATION: Trial no. ' num2str( trial ) ] )
        continue
        
    end
    
    %  samples(:,1)     timestamps of the recording in miliseconds
    %  samples(:,2)     horizontal position of the left eye in degrees
    %  samples(:,3)     vertical position of the left eye in degrees
    %  samples(:,4)     horizontal position of the right eye in degrees
    %  samples(:,5)     vertical position of the right eye in degrees
    
    blinks = zeros( size( samples, 1), 1 );
    
    recording = ClusterDetection.EyeMovRecording.Create(folder, session, ...
        samples, blinks, samplerate);
    [saccades stats] = recording.FindSaccades();
    
    if plotms == 1
        enum = ClusterDetection.SaccadeDetector.GetEnum;
        figure
        subplot(2,2,1)
        plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
        set(gca,'xlim',[0 1],'ylim',[0 100]);
        xlabel('Saccade amplitude (deg)');
        ylabel('Saccade peak velocity (deg/s)');
        
        
        subplot(2,2,[3:4])
        plot(samples(:,1), samples(:,2:end));
        hold
        yl = get(gca,'ylim');
        u1= zeros(size(samples(:,1)))+yl(1);
        u2= zeros(size(samples(:,1)))+yl(1);
        u1((saccades(:,enum.startIndex))) = yl(2);
        u2(saccades(:,enum.endIndex)) = yl(2);
        u = cumsum(u1)-cumsum(u2);
        plot(samples(:,1), u,'k')
        
        xlabel('Time (ms)');
        ylabel('Eye Position (deg)');
        
        legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})
    end
    
    mstrials( 1, trial ) = ~isempty( saccades );
    
    if ~isempty( saccades )
        
        mstrials( 2, trial ) = saccades( 1, 1 ); %time in ms of first microsaccade post fixation
        
        t1 = bhv.CodeTimes{ 1, trial }( find( ismember( bhv.CodeNumbers{ 1, trial }, evts ) ) );
        if ~isempty( t1 ) && any( t1 <= saccades( 1, 1 ) )
            
            mstrials( 3, trial ) = bhv.CodeNumbers{ 1, trial }( bhv.CodeTimes{ 1, trial } == ...
                t1( max( find( t1 <= saccades( 1, 1 ) ) ) ) );
            
        else
            
            mstrials( 3, trial ) = NaN;
            
        end
        
    else
        
        mstrials( 3, trial ) = evts( 1 );
        
    end
    
    clear samples blinks saccades stats recording
    disp( ['Detection Complete: Trial no. ' num2str( trial ) ] )
    
end
end