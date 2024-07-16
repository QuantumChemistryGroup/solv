import std.stdio;
import std.string;
import std.conv;
import std.exception : assertThrown;
import std.typecons;
import std.math;


double mass (string k) {
       double [string] map = [
       "h": 1.0078250,
       "he": 4.0026033,
       "li": 7.0160045,
       "be": 9.0121825,
       "b": 11.0093053,
       "c": 12.0000000,
       "n": 14.0030740,
       "o": 15.9949146,
       "f": 18.9984033,
       "ne": 19.9924391,
       "na": 22.9897697,
       "mg": 23.9850450,
       "al": 26.9815413,
       "si": 27.9769284,
       "p": 30.9737634,
       "s": 31.9720718,
       "cl": 34.9688527,
       "ar": 39.9623831,
       "k": 38.9637079,
       "ca": 39.9625907,
       "sc": 44.9559136,
       "ti": 47.9479467,
       "v": 50.9439625,
       "cr": 51.9405097,
       "mn": 54.9380463,
       "fe": 55.9349393,
       "co": 58.9331978,
       "ni": 57.9353471,
       "cu": 57.9353471,
       "zn": 63.9291454,
       "ga": 68.9255809,
       "ge": 73.9211788,
       "as": 74.9215955,
       "se": 79.9165205,
       "br": 78.9183361,
       "kr": 83.9115064,
       "rb": 84.9117000,
       "sr": 87.9056000,
       "y": 88.9054000,
       "zr": 89.9043000,
       "nb": 92.9060000,
       "mo": 97.9055000,
       "tc": 98.9063000,
       "ru": 101.9037000,
       "rh": 102.9048000,
       "pd": 105.9032000,
       "ag": 106.9050900,
       "cd": 113.9036000,
       "in": 114.9041000,
       "sn": 117.9018000,
       "sb": 120.9038000,
       "te": 129.9067000,
       "i": 126.9004000,
       "xe": 131.9042000,
       "cs": 132.9054290,
       "ba": 137.9050000,
       "la": 138.9061000,
       "hf": 179.9468000,
       "ta": 180.9480000,
       "w": 183.9510000,
       "re": 186.9560000,
       "os": 189.9586000,
       "ir": 192.9633000,
       "pt": 194.9648000,
       "au": 196.9666000,
       "hg": 201.9706000,
       "tl": 204.9745000,
       "pb": 207.9766000,
       "bi": 208.9804000,
       "po": 208.9825000,
       "at": 210.9875000,
       "rn": 222.0175000,
       "fr": 223.0198000,
       "ra": 226.0254000,
       "ac": 227.0278000,
];
  return map [k];
}
