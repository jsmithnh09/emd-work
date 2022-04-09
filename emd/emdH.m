function Y = emdH(IMF)
	% EMD-H.m
	%
	% Original Null Hypothesis EMD-Denoising scheme. Requires the shapiro-wilk
	% testing script. Tests GGD shapes from 1 to 3 in dX = 0.01.
	%
	% [1] Al-Badrawi, Mahdi H., et. al. "A de-noising scheme based on the null
	% hyothesis of intrinsic mode functions." IEEE Signal Processing Letters 23.7
	% (2016): 924 - 928.

	[r, c] = size(IMF);
	if(r > c)
		error('EMDH::incorrect row-transposition of IMFs.');
	end

	% starting parameters.
	alpha = 0.01;
	tail = 0;
	dB = 0.01;
	Limit = 3;
	H(1) = 1;

	for i = 1:r
		[G, pValue, W] = swtest(IMF(i,:), alpha, tail);
		beta = 1;
		R = IMF(i,:);

		while(G == 1) && (beta < Limit)
			sigma = std(R);
			mu = mean(R);
			rho = sqrt( gamma(1/beta)/gamma(3/beta) )*sigma;
			cdf = 0.5 + sign(R - mu).*(gammainc( (abs(R - mu)/rho).^beta, (1/beta) ) );
			beta = beta + dB;
			[G, pValue, W] = swtest(erfinv(2*cdf - 1), alpha, tail);
		end

		H(i) = G;
		shape(i) = beta;
	end

	% Reference selection.
	for d = 2:r
		if H(d) == 1
			m = d;
			break;
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
end  	% EOF
