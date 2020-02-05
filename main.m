clear all;close all;warning off;clc;
addpath([pwd '\Algorithms']);
% run single objective evolutionary algorithms
fobj={  'f10barSI';
        'f25barSI';
        'f37barSI';
        'f60barSI';
        'f72barSI';
        'f120barSI';
        'f200barSI';
        'f942barSI'};    

algo={  'SHAMODE';
        'SHAMODE_WO';};

nvar=[10,8,15,25,16,7,29,59];                % no. of design variables of the test problems
narchive=100;
nbit=5;


for k=1:8 %%%% benchmark
    
    nloop=100;
    nsol=100;
    
    funj=char(fobj(k,:));
    nvari=nvar(k);
    a=zeros(nvari,1);
    b=ones(nvari,1);
       
    for i=1:2  %%%% algorithm
        filename=['rst_' num2str(k,'%03.f') '_' num2str(i,'%03.f') '.mat'];
        for j=1:5     %%%%   optimization run
            RR=[k i j]
            rst(j)=feval(char(algo(i,:)),funj,filename,nloop,nsol,nvari,nbit,narchive,a,b);
        end
        save(filename,'rst','-v7.3');
    end
end   