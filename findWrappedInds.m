function neighboringAngleInds=findWrappedInds(inputInd, MaxInputInd, num_nb)
% find neighboring indices for the three angle parameters, given the number of indices
selIndRange = (inputInd-num_nb):(inputInd+num_nb);
neighboringAngleInds = mod(selIndRange-1, MaxInputInd) + 1;