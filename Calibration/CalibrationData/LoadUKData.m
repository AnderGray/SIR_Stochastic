%%%
%   Loads the Uk data. Cases are daily from 31/01/2020 to 30/03/2020
%   Source: https://www.arcgis.com/home/item.html?id=e5fd11150d274bebaaf8fe2a7a2bda11

Population = 5.5 * 10^6; % About

newCases = ...
[          2
      0
      0
      0
      0
      0
      1
      0
      0
      1
      4
      0
      0
      1
      0
      0
      0
      0
      0
      0
      0
      0
      0
      0
      4
      0
      0
      0
      6
      4
     12
      5
     11
     34
     29
     46
     46
     65
     50
     52
     83
    139
    207
    264
    330
    152
    407
    676
    643
    714
  1035
    665
    967
  1427
  1452
  2129
  2885
  2546
  2433
  2619];

TotalDailyRecovered = ...
[
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
    52
    65
    65
    79
    79
    93
   135
   135
   135
   135
   135
   135
   135
   135
   135
];


TotalDays = length(newCases);

CumInfected = cumsum(newCases);

TotalDailyInfected = CumInfected - TotalDailyRecovered;

