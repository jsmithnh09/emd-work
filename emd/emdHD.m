function Y = emdHD(IMF,X)
	% EMD-Hausdorff Filter. Criterion for selection of reference mode
	% based on the Hausdorff distance between the input signal (prior to decomposition)
	% and the corresponding IMFs. The first IMF with a minimum value after a local maxima
	% is IMFr.
	%
	% SEE ALSO: HausdorffDist.
	%
	% [1] Komaty, Ali, Abdel Boudraa, and Delphine Dare. "EMD-based filtering
	% using the Hausdorff Distance." In Signal processing and Information Technology 
	% (ISSPIT), 2012 IEEE International Symposium on, pp.292-297. IEEE, 2012.

	[r c] = size(IMF);
	if(r > c)
		error('emdHD::incorrect IMF row-transposition.');
	end

	for i = 1:r
		HD(i) = HausdorffDist(X,IMF(i,:));
	end

	% finding local Hausdorff Min. IMF.
	TF = islocalmin(HD);

	% Reference Selection.
	for i = 1:length(TF)
		if TF(i) == 1
			m = i;
		else
			m = 2;
		end
	end

	% Partial Reconstruction.
	rec = 0;
	for j = m:r
		rec = rec + IMF(j,:);
	end
	Y = rec;
end