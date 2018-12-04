function [gwaves, bwaves] = wavecull(waves,wthr)
%
% scan matrix of waveforms and flag outliers
%
%IN: waves = matrix of waveform
%    wthr = similarlity threshold for "good"
%
%OUT: gwaves = indices of wavesforms with similarity to mean >= wthr
%     bwaves = indices of wavesforms with similarity to mean < wthr
%
% AL, janelia, 9.10
%
%[gwaves, bwaves] = wavecull(waves,wthr)
%
%

wmu = mean(waves);
wmu = wmu./norm(wmu);
wlen = size(waves,1);
wavec = zeros(1,wlen);
for k = 1:wlen
    wavek = double(waves(k,:));
    wavec(k) = wmu*wavek'./norm(wavek);
end;
gwaves = find(wavec >= wthr);
bwaves = find(wavec < wthr);

