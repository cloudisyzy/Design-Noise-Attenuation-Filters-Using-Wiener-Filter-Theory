function [x, tempo] = generate_melody(melody_number, Fs)

%
% [x, tempo] = generate_melody(melody_number, Fs)
%	
%	melody_number	- Number between 1 and 5
%	Fs		- Sampling frequency
%	
% 	x		- Generated signal
% 	tempo		- Number of beats per minute
%	
%
%  generate_melody: Generate simple melodies with sinus waves.
%     
%     
%     Mats Bengtsson,  9/2 1998
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch melody_number
  case 1,
    melody=[
      4    -7
      4    -2
      128  -100
      8    -2
      8     0
      4     2
      4    -2
      2     5
      2     2];
    tempo=120;
    
  case 2,
    melody=[
      8     3
      8     0
      8     3
      8     0
      4    -4
      8     3
      8     0
      8     3
      8     0
      4    -4
      8     3
      4     1
      128  -100
      8     1
      4    -2
      128  -100
      8    -2
      2    -4];
    tempo=160;
    
  case 3,
    melody=[
      8  0 
      128  -100
      8  0 
      4  7
      8  7
      16 5
      16 7
      8  3
      16 2
      16 0
      16 2 
      16 3
      4  5
      4  5
      8  7
      8  2
      16 0
      16 -2
      8  -5
      128  -100
      8  -5
      4  0
      8  0
      8  2
      8  0
      16 -2
      16 -3
      16 -2
      16 -3
      2  -5];
    tempo=80;
    
  case 4,
    melody=[
      2  -2
      4  2
      2  -2
      4  -5
      2  -9
      1  0
      4  -4
      2  -3
      4  0
      2  -3
      4  -7
      2  -10
      4  -2];
    tempo=360;
    
  case 5,
    melody=[
      4  -2
      8  -2
      8  -7
      4  -2
      8  -2
      8  -7
      8  -2
      8  -7 
      8  -2
      8   2
      2   5
      4   3
      8   3
      8   0
      4   3
      8   3
      8   0
      8   3
      8   0
      8  -3
      8   0
      2  -7
      ];
    tempo=120;
      
end

x=sequencer(melody,tempo,Fs);
