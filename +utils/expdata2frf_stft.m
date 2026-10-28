function G = expdata2frf_stft(matfile, varargin)
%EXPDATA2FRF  FRF of a single-sine sweep measurement (ridge / sweep-synchronous estimate).
%
%   G = utils.expdata2frf('file.mat')
%   G = utils.expdata2frf('file.mat','InputVar','u','OutputVar','y','Fs',1000)
%
%   Instead of Y(w)/U(w) on the whole record (which mixes in transients and
%   distortion occurring at other times), the FRF is estimated at the instant
%   the sweep passes each frequency, with a sliding windowed single-frequency
%   DFT of input and output.
%
%   Returns an frd object (frequency in rad/s), ny outputs x 1 input.
%
%   Options (name-value):
%     'InputVar'        name of the input variable in the .mat   (auto-detected)
%     'OutputVar'       name of the output variable(s)            (auto-detected)
%     'Fs'              sampling frequency in Hz                  (from file, else 1000)
%     'Window'          window length in samples                  (1024)
%     'Hop'             hop between estimates in samples          (64)
%     'FreqRange'       [fmin fmax] in Hz kept in the result      ([10 400])
%     'ActiveThreshold' fraction of max |u| defining sweep on/off (0.01)
%     'Plot'            true -> bode plot                         (false)
%
%   Assumes a LINEAR sweep (as in your spectrogram). Output y may be a matrix
%   (one column per accelerometer).

p = inputParser;
addParameter(p,'InputVar','');
addParameter(p,'OutputVar','');
addParameter(p,'Fs',[]);
addParameter(p,'Window',1024);
addParameter(p,'Hop',64);
addParameter(p,'FreqRange',[10 500]);
addParameter(p,'ActiveThreshold',0.01);
addParameter(p,'Plot',false);
addParameter(p,'itmax',[]);
addParameter(p,'itmin',1);
parse(p,varargin{:});
o = p.Results;

%% Load data
[u,y]=utils.openIOtemp(matfile,o.itmin,o.itmax);
if isvector(y), y = y(:); end
if size(y,1) ~= numel(u) && size(y,2) == numel(u), y = y.'; end
assert(size(y,1) == numel(u), 'Input and output have different lengths.');

%% Sampling frequency
fs = o.Fs;
if isempty(fs)
    if     isfield(S,'fs'), fs = S.fs;
    elseif isfield(S,'Fs'), fs = S.Fs;
    elseif isfield(S,'Ts'), fs = 1/S.Ts;
    elseif isfield(S,'dt'), fs = 1/S.dt;
    elseif isfield(S,'t') && numel(S.t) > 1, fs = 1/mean(diff(S.t));
    elseif isfield(S,'time') && numel(S.time) > 1, fs = 1/mean(diff(S.time));
    else
        fs = 1000;
        warning('expdata2frf:fs','No sampling info in file, assuming Fs = 1000 Hz.');
    end
end

N = o.Window;  hop = o.Hop;
nS = numel(u);
w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);          % Hann window

%% Active part of the sweep
env = movmax(abs(u), 200);
act = find(env > o.ActiveThreshold*max(abs(u)));
i0 = act(1);  i1 = act(end);

%% Sweep law f(t) = a + b*t from the ridge of the input spectrogram
[Su, f, t] = spectrogram(u, hann(N), N-hop, N, fs);
[~, k] = max(abs(Su), [], 1);
fr = f(k).';  t = t(:).';
inAct = t >= (i0 + N/2)/fs & t <= (i1 - N/2)/fs;
fr_a = fr(inAct);  t_a = t(inAct);
span = max(fr_a) - min(fr_a);
core = fr_a > min(fr_a) + 0.05*span & fr_a < max(fr_a) - 0.05*span;
c = polyfit(t_a(core), fr_a(core), 1);       % c(1) = Hz/s, c(2) = Hz at t = 0

%% FRF at the instant the sweep passes each frequency
ic = (N/2+1) : hop : (nS - N/2);             % window centres (samples)
tc = (ic-1)/fs;
fc = polyval(c, tc);
keep = fc >= o.FreqRange(1) & fc <= o.FreqRange(2) ...
     & ic >= i0 + N/2 & ic <= i1 - N/2;
ic = ic(keep);  fc = fc(keep);  K = numel(ic);

tl  = (-N/2 : N/2-1)'/fs;
idx = ic + (-N/2 : N/2-1)';                  % N x K sample indices
E   = exp(-1j*2*pi*tl*fc) .* w;              % N x K

Uc = sum(u(idx) .* E, 1).';                  % K x 1
ny = size(y,2);
H  = zeros(K, ny);
for m = 1:ny
    ym = y(:,m);
    H(:,m) = sum(ym(idx) .* E, 1).' ./ Uc;
end

%% Build frd
resp = reshape(H.', ny, 1, K);
G = frd(resp, 2*pi*fc,0.001);
% [~, stem] = fileparts(matfile);
% G.Name = stem;
% G.InputName = 'u';
% if ny == 1, G.OutputName = 'y'; else, G.OutputName = compose('y%d', 1:ny)'; end

if o.Plot
    bode(G); grid on
end
end

%% ---------- helpers ----------
function S = flattenStruct(S)
% If the .mat holds one struct (e.g. 'data'), use its fields.
fn = fieldnames(S);
if numel(fn) == 1 && isstruct(S.(fn{1}))
    S = S.(fn{1});
end
end

function name = pick(S, requested, candidates, label)
fn = fieldnames(S);
if ~isempty(requested)
    assert(isfield(S, requested), 'Variable "%s" not found. Available: %s', ...
        requested, strjoin(fn', ', '));
    name = requested;
    return
end
for i = 1:numel(candidates)
    if isfield(S, candidates{i}), name = candidates{i}; return, end
end
error('expdata2frf:var', ['Could not auto-detect the %s variable. Available: %s\n' ...
    'Pass it with ''InputVar'' / ''OutputVar''.'], label, strjoin(fn', ', '));
end
