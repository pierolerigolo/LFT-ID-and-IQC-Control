function [A,B,C,D] = Gdm(dm,c,k,Ts,ps)
%GDM Grey-box model for utils.Gdm_greyestid (idgrey ODE function)
%
%   1 x n system, channel i = s^2/(dm(i) s^2 + c s + k), Tustin-discretised
%   (must match utils.Gdm_greyestid). Each channel keeps 2 states, with a
%   fixed sign convention per channel and the phase sign ps applied.
%
    n  = numel(dm);
    Gi = cell(n,1);
    for i = 1:n
        Gi{i} = ss(tf([1 0 0],[dm(i) c k]));    % always 2 states
    end
    G = horzcat(Gi{:});                      % 1 x n, 2n states, no reduction
    if Ts ~= 0
        G = c2d(G,Ts,'tustin');              % keeps the order; makes the mass LFT exact
    end
    [A,B,C,D] = ssdata(G);

    % sign convention per input channel j (states 2j-1:2j)
    for j = 1:n
        idx = 2*j-1:2*j;
        cj  = 0;
        nz  = find(C(1,idx) ~= 0, 1);
        if ~isempty(nz), cj = C(1,idx(nz)); end
        if ~(D(1,j) < 0 || (D(1,j) == 0 && cj < 0))
            B(idx,j) = -B(idx,j);
            D(1,j)   = -D(1,j);
        end
    end

    % phase convention, fixed extra argument (scalar or 1 x n)
    B = B.*ps;
    D = D.*ps;
end
