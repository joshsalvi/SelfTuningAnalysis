classdef DataImport
    properties
        Fs
        pre
        pulse
        cyclesdur
        cyclesdurlength
        post
        logdata
        comments
        ntrace
        numavg
        raw
        nonraw
        Xd
        Xo
        Fe
        Xd_pulse
        Xo_pulse
        Fe_pulse
        Xd_split
        Xo_split
        kv
        datapath
        tvec
    end
    methods
        function savedata(obj,datapath)
            disp('Saving...');
            savefile = sprintf('%s%s',datapath,'Extracted Data-classfile.mat');
            save(savefile,'obj');
            disp('Finished.');
        end
        function obj = initialization(obj,datapath,Fs,pre,pulse,cyclesdur,cyclesdurlength,post,ntrace,numavg,raw,nonraw)
            obj.datapath = datapath;
            obj.Fs = Fs;
            obj.pre = pre;
            obj.pulse = pulse;
            obj.cyclesdur = cyclesdur;
            obj.cyclesdurlength = cyclesdurlength;
            obj.post = post;
            obj.ntrace = ntrace;
            obj.numavg = numavg;
            obj.raw = raw;
            obj.nonraw = nonraw;
        end
        function obj = logfileimport(obj,logdata)
            obj.logdata = logdata;
            if isempty('logdata.textdata(isnan(logdata.data(:,8))==0,3)')==0
                obj.comments = logdata.textdata(isnan(logdata.data(:,8))==0,3); % import comments
            else
                obj.comments = 'no comments';
            end
        end
        function obj = timeseries(obj,tvec,Xd,Xo,Fe,Xd_pulse,Xo_pulse,Fe_pulse,Xd_split,Xo_split,kv)
            obj.tvec = tvec;
            obj.Xd = Xd;
            obj.Xo = Xo;
            obj.Fe = Fe;
            obj.Xd_pulse = Xd_pulse;
            obj.Xo_pulse = Xo_pulse;
            obj.Fe_pulse = Fe_pulse;
            obj.Xd_pulse = Xd_split;
            obj.Xo_split = Xo_split;
            obj.kv = kv;
        end
                     
    end
end
