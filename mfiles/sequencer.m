function x = sequencer(descr, tempo, Fs)

%
% x = sequencer(descr, tempo, Fs)
%	
%	descr	- Matrix where each row describes a tone
%		  in the form [length, pitch]
%		length	- Note duration as fraction of a wholenote.
%		pitch	- Pitch in semitones over mid A
%	tempo	- Number of fourth notes (quavers)/minute
%	Fs	- Sampling frequency
%	
% 	x	- Sample music.
%	
%
%  sequencer:
%     Toy to generate sampled music from a description.
%     
%     Mats Bengtsson,  8/1 1998
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note durations in number of samples:
durations = round(60/tempo*4*Fs./descr(:,1));
N=sum(durations); % Number of samples

% Normalized frequency of each tone:
% (a=440 Hz, c is 9 semitones lower, 12 semitones form one octave)
frequencies= 440/Fs*pow2((descr(:,2))/12);

% Instantanous frequency at each sample:
q=zeros(N,1);
t=0;
for n=1:size(descr,1);
  q(t+1:t+durations(n))=frequencies(n)*ones(durations(n),1);
  t=t+durations(n);
end

x=sin(2*pi*cumsum(q));
