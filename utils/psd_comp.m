function psd_comp(a,b,c,d)

    [pa, f] = pwelch(a, hamming(64), 32, 2048, 1, 'twosided');
    N = length(f);
    halfN = floor(N/2);
    f = f(1:halfN);
    pa = pa(1:halfN);

    [pb, ~] = pwelch(b, hamming(64), 32, 2048, 1, 'twosided');
    pb = pb(1:halfN);

    [pc, ~] = pwelch(c, hamming(64), 32, 2048, 1, 'twosided');
    pc = pc(1:halfN);

    [pd, ~] = pwelch(d, hamming(64), 32, 2048, 1, 'twosided');
    pd = pd(1:halfN);

    semilogy(f, pa, 'b', 'LineWidth',1.5); hold on 
    semilogy(f, pb, ':k', 'LineWidth',1.5); hold on 
    semilogy(f, pc, '-.r', 'LineWidth',1.5); hold on 
    semilogy(f, pd, '--g', 'LineWidth',1.5); 

end

