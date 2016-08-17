function out = amtafcexp(flag,par,varargin)
%AMTAFTEXP Runs experimental trials for binaural singal detection with a model
% 
%   Usage:   out = amtafcexp(flag,varargin) 
% 
%   Input parameters:
%       flag:       specifies the chosen execution 
%       par:        already set experimental parameters
%                   
%       The following flags are possible:
%           'expinit':      sets all parameters for the experimental run
%           'modelinit':    sets all model parameters
%           'signalinit':   sets all signal parameters
%           'decisioninit': sets all decision parameters
%           'run':          runs the experiment and lets the model decide
%
%   run 'result = amtafcexp('run',par,'plot')' to generate experimental plot


switch flag
    case 'expinit'
        definput.keyvals.intnum = [];
        definput.keyvals.rule = [];
        definput.keyvals.expvarstepstart = [];
        definput.keyvals.expvarsteprule = [];
        definput.keyvals.stepmin = [];
        definput.keyvals.expvarstart = [];
        
        [~,kvexp]=ltfatarghelper({},definput,varargin{:});
             
       % out = setstructfields(kvexp, par);
       out = par;
       out.exp = kvexp;
        
    case 'modelinit'
        definput.keyvals.name=[];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];
        definput.keyvals.outputs = [];
        
        [~,kvmodel]=ltfatarghelper({},definput,varargin{:});
        
        % check what outputs of model are needed
        outputnumber = nargout(kvmodel.name);
        
        callmodelstring = [];

        for outputcounter = 1:outputnumber
            if any(kvmodel.outputs==outputcounter)
                modelstring = sprintf('modelout.par%i',outputcounter);
                if isempty(callmodelstring)
                    callmodelstring = ['[' modelstring '{interval_num}'];
                else
                    callmodelstring =  [callmodelstring ',' modelstring '{interval_num}'];
                end
            else
                callmodelstring = [callmodelstring ',~'];           
            end
        end
        
        modelinputs = struct2cell(kvmodel);
                
        % delete name of model function & outputs
        modelinputs = modelinputs(2:end-1);

        % delete empty cells
        modelinputs = modelinputs(~cellfun('isempty',modelinputs));
        
        callmodelstring = [callmodelstring ']=' kvmodel.name '('];
        
        for modelinputscounter = 1:length(modelinputs)
            callmodelstring = [callmodelstring num2str(modelinputs{modelinputscounter}) ','];
        end
        callmodelstring(end:end+1) = ');';
        
        out = par;
        out.model = kvmodel;
        out.callstrings.model = callmodelstring;
        
        %TODO CALLMODELSTRING
        
    case 'signalinit'
        definput.keyvals.name= [];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];
        definput.keyvals.input11 = [];
        definput.keyvals.input12 = [];
        definput.keyvals.input13 = [];
        definput.keyvals.input14 = [];
        definput.keyvals.input15 = [];

        
        [~,kvsignal]=ltfatarghelper({},definput,varargin{:});
        
        out = par;
        out.signal = kvsignal;
    
    case 'decisioninit'
        definput.keyvals.name= [];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];
        
        definput.flags.plot = {'noplot','plot'};
        
        [~,kvdecision]=ltfatarghelper({},definput,varargin{:});
        
        out = par;
        out.decision = kvdecision;
        
    case 'run'
        
        definput.flags.plot = {'noplot','plot'};
        
        [flags,~]=ltfatarghelper({},definput,varargin);
        
        % find experimental variable and inttyp variable
        sigparnames = fieldnames(par.signal);
        for count = 1: length(sigparnames)
            name = sigparnames(count);
            if strcmp(getfield(par.signal, name{:}),'expvar')
                experimentvar = name{:};
            elseif strcmp(getfield(par.signal, name{:}),'inttyp')
                inttypvar = name{:};
            end
        end
        
        stepsize = par.exp.expvarstepstart;
        truecounter = 0;
        par.signal.(experimentvar) = par.exp.expvarstart;          
        expparvalue = [];
        downturn = 0;
        upturn = 1;
        turncounter = 0;
        lastturn = [];
        checkmodelout = 0;
        condition = 1;
        
        while condition
            
            for interval_num=1:par.exp.intnum
                
                if interval_num == 1
                    par.signal.(inttypvar) = 'target'; 
                else
                    par.signal.(inttypvar) = 'reference';
                end
                          
                signalinputs = struct2cell(par.signal);
                
                % delete name of signal function
                signalinputs = signalinputs(2:end);
                
                % delete empty cells
                signalinputs = signalinputs(~cellfun('isempty',signalinputs));
                
                % call signalfunction
                testsignal = feval(par.signal.name,signalinputs{:});
                
                % find experimental signal variable
                modelparnames = fieldnames(par.model);
                for count = 1: length(modelparnames)
                    name = modelparnames(count);
                    if strcmp(getfield(par.model, name{:}),'expsignal')
                        experimentsignal = name{:};
                        break
                    end
                end
                par.model.(experimentsignal) = testsignal;                           
                
                % call model
                par.callstrings.model = strrep(par.callstrings.model,'expsignal', 'testsignal');
                eval(par.callstrings.model);
                
            end
                
            % find model output variable, only at the first time
            if checkmodelout == 0
                modeloutvarcount = 1;
                decisionparnames = fieldnames(par.decision);
                for count = 1:length(decisionparnames)
                    name = decisionparnames(count);
                    if strcmp(getfield(par.decision, name{:}),'modelout')
                        modeloutvar{modeloutvarcount} = name{:};
                        modeloutvarcount = modeloutvarcount+1;
                    end
                end
                checkmodelout = 1;
            end

            %set model outputs to decsion inputs
            for count = 1:modeloutvarcount-1
                modeloutname = sprintf('par%i',count);
                par.decision.(modeloutvar{count}) = modelout.(modeloutname);
            end

            decisioninputs = struct2cell(par.decision);

            % delete name of model function & outputs
            decisioninputs = decisioninputs(2:end);

            % delete empty cells
            decisioninputs = decisioninputs(~cellfun('isempty',decisioninputs));

            % call decision
            decision = feval(par.decision.name,decisioninputs{:});
            
            % store expparvalue
            expparvalue = [expparvalue par.signal.(experimentvar)];
            
            % count reversals
            % wrong answers are par.exp.rule(2)and no low point reversal
            if decision ~= 1 && wrongcounter == par.exp.rule(2)-1 && downturn == 0
                turncounter = turncounter + 1;
                downturn = 1;
                upturn = 0;
                if stepsize == par.exp.stepmin(1) 
                    lastturn = [lastturn par.signal.(experimentvar)];
                end

             % right answers are par.exp.rule(1) 
             % and no high point reversal
            elseif decision == 1 && truecounter == par.exp.rule(1)-1 && upturn == 0
                turncounter = turncounter + 1;
                downturn = 0;
                upturn = 1;
                if stepsize == par.exp.stepmin(1)
                    lastturn = [lastturn par.signal.(experimentvar)];
                end
            end
       
            % change stepsize after par.exp.expvarsteprule(2) reversals, if stepsize
            % is not already par.exp.stepmin(1) dB
            if turncounter == par.exp.expvarsteprule(2) && truecounter == 1 && ...
                stepsize ~= par.exp.stepmin(1)
                stepsize = stepsize * par.exp.expvarsteprule(1);
                turncounter = 0;
            end
            
            
            if decision == 1
                wrongcounter = 0;
                truecounter = truecounter + 1;                
            else
                truecounter = 0;
                wrongcounter = wrongcounter +1;
            end
            
            % amount of right answers is par.exp.rule(1)
            if truecounter == par.exp.rule(1)
                par.signal.(experimentvar) = par.signal.(experimentvar) - stepsize;
                truecounter = 0;
            end 
            
            % amount of wrong answers is par.exp.rule(2)
            if wrongcounter == par.exp.rule(2)
                par.signal.(experimentvar) = par.signal.(experimentvar) + stepsize;
                wrongcounter = 0;
            end 
            
            % amount of reversals is par.stepmin(2)
            if size(lastturn,2) == par.exp.stepmin(2)
                threshold = median(lastturn);
                threshstd = std(lastturn,1);
                condition = 0;
                out = [threshold threshstd expparvalue];
            end
   
        end
        %clear persistent variables
        clear (par.decision.name);
        
        if flags.do_plot
            figure
            for plotcounter = 1:(size(expparvalue,2)-1)
                if expparvalue(plotcounter) < expparvalue(plotcounter+1)
                    style = 'or';
                    stem(plotcounter,(expparvalue(plotcounter)),'or',...
                        'LineStyle','none')
                    hold on
                else
                    style = '+g';
                    stem(plotcounter,(expparvalue(plotcounter)),'+g',...
                        'LineStyle','none')
                    hold on
                end
            end
            
            % special case last entry 
            if decision == 1
                stem(plotcounter+1,(expparvalue(plotcounter+1)),'+g',...
                        'LineStyle','none')
            else
                stem(plotcounter+1,(expparvalue(plotcounter+1)),'or',...
                        'LineStyle','none')
            end 
            title(['Threshold (Median): ' num2str(threshold) ...
                'dB, Std: ' num2str(threshstd) 'dB'])
            
        end
        
end
