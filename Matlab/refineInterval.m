function tt = refineInterval(Interval, numRefine, Ngauss)
    

    toBeAdded = zeros(2*length(Interval)-2, numRefine);
    
    toBeAdded(1, 1) = (Interval(1) + Interval(2))/2;
    for j = 2:numRefine
       toBeAdded(1, j) = (Interval(1) +  toBeAdded(1, j-1))/2;
    end
    
    
    for i = 2:(length(Interval)-1)
       
        toBeAdded(i, 1) = (Interval(i-1) + Interval(i))/2;
        for j = 2:numRefine
            toBeAdded(i, j) = (Interval(i) + toBeAdded(i, j-1))/2;
        end
        toBeAdded(i + length(Interval) - 2, 1) = (Interval(i) + Interval(i+1))/2;
        for j = 2:numRefine
            toBeAdded(i + length(Interval) - 2, j) = (Interval(i) + toBeAdded(i + length(Interval) - 2, j-1))/2;
        end
    end

    
    toBeAdded(end, 1) = (Interval(end - 1) + Interval(end))/2;
    
     for j = 2:numRefine
       toBeAdded(end, j) = (Interval(end) +  toBeAdded(end, j-1))/2;
    end

        
    tt = sort(unique([reshape(toBeAdded, 1, []), Interval]));  
    
       
    %[tt, ww] = makeGaussianNodes(t, Ngauss);

end




function [tNodes, wNodes] = makeGaussianNodes(list, numNodes)
   Nnodes = length(list) - 1;
   nodesAtList = zeros(Nnodes);
   if Nnodes == length(numNodes)
       nodesAtList = numNodes;
   else
      nodesAtList(:) = numNodes(1,1); 
   end
    [tNodes, wNodes] = lgwt(nodesAtList(1), list(1), list(2));
   for i = 2:Nnodes
        [a, b] = lgwt(nodesAtList(i), list(i), list(i + 1));
        tNodes = [tNodes; a];
        wNodes = [wNodes; b];
   end
end
