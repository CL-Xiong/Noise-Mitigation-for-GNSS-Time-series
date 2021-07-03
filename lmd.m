function [pf,a,s] = lmd(x)

c = x(:)'; % copy of the input signal (as a row vector)
N = length(x);
range= max(x)-min(x); %ÐÅºÅµÄ·åÖµ£¬
%-------------------------------------------------------------------------
% loop to decompose the input signal into successive PF

pf = []; % Matrix which will contain the successive IMF, and the residue
a = [];
s = [];
count = 0; %
h = c; % at the beginning of the sifting process, h is the signal
while (1) 
   if (max(h)-min(h))<range*0.001   % if the amplitude is much less than the amplitude of original signal ,than return
       pf=[pf,h']; % save the residue
       return
   end
   an=ones(N,1);   % save the envelopes an=an1*an2*an3*...*an
   while (1) % loop to find a PF      
      % find local max/min points
      d = diff(h); % approximate derivative
      maxmin = []; % to store the optima (min and max without distinction so far)
      for i = 1:N-2
         if (d(i)==0) && (i~=1) && (sign(d(i-1))~=sign(d(i+1)))                % we are on a zero
            maxmin = [maxmin, i];
         elseif (d(i)~=0)&&(sign(d(i))~=sign(d(i+1)))   % we are straddling a zero so
            maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
         end
      end
      
      if size(maxmin,2) < 2 % then it is the residue
          pf = [pf,h']; % save the residue
          return 
      end
      maxmin = [1,maxmin,N];
      mi = [];
      ai = [];
      for i = 1:length(maxmin)-1
          m_seg = ones(1,maxmin(i+1)-maxmin(i));   % calculate a segment of mean of envelope and magnitude. 
          a_seg = ones(1,maxmin(i+1)-maxmin(i));
          m_seg = m_seg*(h(maxmin(i))+h(maxmin(i+1)))/2;
          a_seg = a_seg*abs((h(maxmin(i))-h(maxmin(i+1))))/2;
          mi = [mi,m_seg];
          ai = [ai,a_seg];
      end  
      mi = [mi,mi(length(mi))];% add a element
      ai = [ai,ai(length(ai))];
      for i = 2:length(maxmin)-1
          mi(maxmin(i)) = (mi(maxmin(i))+mi(maxmin(i)-1))/2;
          ai(maxmin(i)) = (ai(maxmin(i))+ai(maxmin(i)-1))/2;
      end
      span = 0; %determine the span of the average smoothing
      for i = 2:length(maxmin)
          if (maxmin(i)-maxmin(i-1))>span
              span = maxmin(i)-maxmin(i-1);
          end
      end
      span = round(span/3); 
      if span<3   % span should be bigger than a number
          span = 3;
      end
      while(1)      % average smoothing until no two successive has the same value;
        mi = smooth(mi,span);
        ai = smooth(ai,span);
        I = find(diff(mi)==0);
        if length(I)==0
            break;
        end
    end
    an = an.*ai;
    h = h-mi';
    si = h./ai';      % frequency modulated signal
    h = si;
    delta = 0.1; % in practical,for stop condition 
    if ( min(ai)>1-delta && max(delta)<1+delta)  % if ai is near 1,the a purely frequency modulate signal is gained.
      for i = 1:length(si)
          if si(i)>1
              si(i) = 1;
          end
          if si(i)<-1
              si(i)=-1;
          end
      end
      pf = [pf,an.*si'];
      a = [a,an];
      s = [s,si'];
      c = c-an'.*si; % the original signal minus the product function
      h = c;
      break;
  end
end  
end
