function T = perftable(Ks, names, GFRF, omega_d, opts)
%PERFTABLE Attenuation |S| at the disturbance frequency on every operating point
%
%   T = utils.perftable(Ks, names, GFRF, omega_d)
%   T = utils.perftable(..., 'realmass', m, 'G_hat', G_hat, 'file', 'plots/x.csv')
%
%   Ks     cell of controllers, names their column labels (string array)
%   GFRF   measured FRFs, one channel / model per operating point
%   One row per operating point (plus the nominal model when G_hat is given),
%   one column per controller: |S(j omega_d)| = |1/(1 + K G_i)| in dB.
%   The table is printed, written to opts.file (csv, if given) and copied to
%   the clipboard as tab-separated text: paste it in Keynote/Numbers/Excel
%   and it becomes a table.
%
arguments
    Ks cell
    names string
    GFRF
    omega_d (1,1) double
    opts.realmass double = []
    opts.G_hat = []
    opts.file string = ""
    opts.digits (1,1) double = 1
    opts.summary (1,1) logical = true   % add worst / best / mean / spread rows (operating points only)
    opts.GFRFs cell = {}                % optional: one GFRF per controller (e.g. with / without a
                                        % compensator), same operating points; default GFRF for all
end
    nOP = size(GFRF, 2);          % one operating point per channel (as in the c scripts)
    Gw  = zeros(nOP, 1);
    for i = 1:nOP
        Gw(i) = squeeze(interp(GFRF(i), omega_d).ResponseData);
    end
    rows = "OP" + string((1:nOP).');
    if ~isempty(opts.realmass)
        rows = rows + " (" + compose("%.2f", opts.realmass(:)) + " kg)";
    end
    if ~isempty(opts.G_hat)
        Gw(end+1) = squeeze(interp(opts.G_hat, omega_d).ResponseData);
        rows(end+1) = "Nominal";
    end

    S = zeros(numel(Gw), numel(Ks));
    for j = 1:numel(Ks)
        Gj = Gw;
        if numel(opts.GFRFs) >= j && ~isempty(opts.GFRFs{j})
            for i = 1:nOP, Gj(i) = squeeze(interp(opts.GFRFs{j}(i), omega_d).ResponseData); end
        end
        S(:, j) = mag2db(abs(1./(1 + squeeze(freqresp(Ks{j}, omega_d))*Gj)));
    end
    if opts.summary
        Sop  = S(1:nOP, :);                          % measured operating points only
        S    = [S; max(Sop, [], 1); min(Sop, [], 1); mean(Sop, 1); max(Sop, [], 1) - min(Sop, [], 1)];
        rows = [rows; "Worst (max dB)"; "Best (min dB)"; "Mean"; "Spread (max - min)"];
    end
    S = round(S, opts.digits);

    head = "|S| at " + compose("%.0f", omega_d/2/pi) + " Hz (dB) - " + names(:).';
    T = array2table(S, 'VariableNames', cellstr(head), 'RowNames', cellstr(rows));
    disp(T)

    if strlength(opts.file) > 0
        writetable(T, opts.file, 'WriteRowNames', true);
    end
    fmt = "%." + opts.digits + "f";
    lines = ["Operating point", head; rows, compose(fmt, S)];
    txt = strjoin(join(lines, char(9), 2), newline);
    try
        clipboard('copy', char(txt));
        fprintf('Table copied to the clipboard (tab separated): paste it in Keynote.\n');
    catch
    end
end
