function Y = cemd(IMF)
	% Conventional EMD (CEMD)
	%
	% EMD filter that measures the power present within the first noise
	% contaminated IMF to formulate a threshold in selecting the reference IMF.
	%
	% [1] P. Flandrin, P. Goncalves, and G. Rilling,
	% "EMD Equivalent Filter Banks, from interpretation to applications."
	% in Hilbert-Huang Tranform: Introduction and Applications,
	% vol. 5, Singapore: World Scientific, 2005.

	[r, c] = size(IMF);
	if(r > c)
		error('CEMD:: IMFs are not row transposed.');
	end

	E1 = sum(IMF(1,:).^2);

	for i = 1:r-1
		E(i) = (E1/0.719)*(2.01)^-(i+1);
		% T95(i) = 2^(0.474*i-2.449)+log2(E(i));
		T99(i) = 2^(0.46*i-1.919)+log2(E(i));
		Ey(i) = log2(sum(IMF(i+1,:).^2));
	end

	% Reference Selection.
	for dc = 1:length(Ey)
		if Ey(dc) > T99(dc)
			mc = dc;
			break;
		else
			mc = 2;
		end
	end

	% Partial reconstruction
	rec = 0;
	for j = mc:r
		rec = rec + IMF(j,:);
	end
	Y = rec;
end