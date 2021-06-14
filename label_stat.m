function [total_stat,chan_stat,TP,FP,FN]=label_stat(test,ref,tol)

% INPUTS:------------------------------------------------------------------
%   test... structure of tested events
%           test.pos (column vector) - position of event in second
%           test.chan (column vector) - channel of event
%           test.con (column vector) - condition of event (not used)
%   ref... structure of reference events
%           ref.pos (column vector) - position of event in second
%           ref.chan (column vector) - channel of event
%           ref.con (column vector) - condition of event (not used)
%   tol... tolerance boundaries between test and reference-event for match
%           in second (recomended 0.1)
%           EXAMPLE: test.pos=10.1 +/- 0.1; 
%                    ref.pos=10.05;
%                    tested event is accepted as true positive (TP) 
%
% OUTPUTS:-----------------------------------------------------------------
%   total_stat... [TP,FP,FN,SEN,SEL] 
%                  TP (num. of true positive) - match
%                  FP (num. of false positive) - false detection
%                  FN (num. of false negative) - unmarked
%   chan_stat... matrix ch x 5, where ch is nuber of channels. 
%                Each row contains stat values [TP,FP,FN,SEN,SEL]  
%   TP... structures of events of true positives events
%           TP.pos (column vector) - position of event in second
%           TP.chan (column vector) - channel of event
%           TP.con (column vector) - condition of event (not used)
%   FP... false postive, same structure as TP
%   FN... false postive, same structure as TP
%
% MADE BY ISARG - Radek Janca, 2013
% http://sami.fel.cvut.cz/isarg/


%%
num_test=length(test.pos);
num_ref=length(ref.pos);

% comparison and sorting to class FP, TP and FN
TP.con=[];
TP.chan=[];
TP.pos=[];

FP.con=[];
FP.chan=[];
FP.pos=[];

while ~isempty(test.pos)    
    idx=find(ref.pos>test.pos(1)-tol & ref.pos<test.pos(1)+tol & ref.chan==test.chan(1)); % comparison test to ref
       
    if isempty(idx) % no match - FP
       FP.con=[FP.con; test.con(1)];
       FP.chan=[FP.chan; test.chan(1)];
       FP.pos=[FP.pos; test.pos(1)];
       
       test.con(1)=[]; % removing from test events
       test.chan(1)=[];
       test.pos(1)=[];
    
    else % match - TP
        
       idx=idx(1);
       TP.con=[TP.con; test.con(1)];
       TP.chan=[TP.chan; test.chan(1)];
       TP.pos=[TP.pos; test.pos(1)];
       
       test.con(1)=[]; % removing from test events
       test.chan(1)=[];
       test.pos(1)=[];
       
       ref.con(idx)=[]; % removing from ref events
       ref.chan(idx)=[];
       ref.pos(idx)=[];
    end 
end

FN=ref; % remainder of ref events is FN

%%
total_stat(1,1)=length(TP.pos);
total_stat(1,2)=length(FP.pos);
total_stat(1,3)=length(FN.pos);
total_stat(1,4)=length(TP.pos)/(length(TP.pos)+length(FN.pos));
total_stat(1,5)=length(TP.pos)/(length(TP.pos)+length(FP.pos));

disp(['num of TP: ' num2str(total_stat(1,1))])
disp(['num of FP: ' num2str(total_stat(1,2))])
disp(['num of FN: ' num2str(total_stat(1,3))])
disp(['---------------'])
disp(['sensitivity SEN= ' num2str(total_stat(1,4),'%0.3f')])
disp(['selectivity SEL= ' num2str(total_stat(1,5),'%0.3f')])

if sum([length(TP.pos) length(FP.pos)])~=num_test
    warning(['TP+FP: ' num2str(sum([total_stat(1,1) total_stat(1,2)])) ' (in test ' num2str(num_test) ')'])
end
if sum([length(TP.pos) length(FN.pos) ])~=num_ref
    warning(['TP+FN: ' num2str(sum([total_stat(1,1) total_stat(1,3)])) ' (in ref ' num2str(num_ref) ')'])
end

%% each channel satatistics - Tp,FP,FN,SEN,SEL of each channel

for ch=1:max([TP.chan;FN.chan;FP.chan])
    chan_stat(ch,1)=sum(TP.chan==ch);
    chan_stat(ch,2)=sum(FP.chan==ch);
    chan_stat(ch,3)=sum(FN.chan==ch);
    chan_stat(ch,4)=sum(TP.chan==ch)/(sum(TP.chan==ch)+sum(FN.chan==ch));
    chan_stat(ch,5)=sum(TP.chan==ch)/(sum(TP.chan==ch)+sum(FP.chan==ch));
end


