$offlisting
*  
*  Equation counts
*      Total        E        G        L        N        X        C        B
*        401      401        0        0        0        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*        507      507        0        0        0        0        0        0
*  FX      7
*  
*  Nonzero counts
*      Total    const       NL      DLL
*       2002      802     1200        0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19
          ,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36
          ,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53
          ,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70
          ,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87
          ,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100,x101,x102,x103
          ,x104,x105,x106,x107,x108,x109,x110,x111,x112,x113,x114,x115,x116
          ,x117,x118,x119,x120,x121,x122,x123,x124,x125,x126,x127,x128,x129
          ,x130,x131,x132,x133,x134,x135,x136,x137,x138,x139,x140,x141,x142
          ,x143,x144,x145,x146,x147,x148,x149,x150,x151,x152,x153,x154,x155
          ,x156,x157,x158,x159,x160,x161,x162,x163,x164,x165,x166,x167,x168
          ,x169,x170,x171,x172,x173,x174,x175,x176,x177,x178,x179,x180,x181
          ,x182,x183,x184,x185,x186,x187,x188,x189,x190,x191,x192,x193,x194
          ,x195,x196,x197,x198,x199,x200,x201,x202,x203,x204,x205,x206,x207
          ,x208,x209,x210,x211,x212,x213,x214,x215,x216,x217,x218,x219,x220
          ,x221,x222,x223,x224,x225,x226,x227,x228,x229,x230,x231,x232,x233
          ,x234,x235,x236,x237,x238,x239,x240,x241,x242,x243,x244,x245,x246
          ,x247,x248,x249,x250,x251,x252,x253,x254,x255,x256,x257,x258,x259
          ,x260,x261,x262,x263,x264,x265,x266,x267,x268,x269,x270,x271,x272
          ,x273,x274,x275,x276,x277,x278,x279,x280,x281,x282,x283,x284,x285
          ,x286,x287,x288,x289,x290,x291,x292,x293,x294,x295,x296,x297,x298
          ,x299,x300,x301,x302,x303,x304,x305,x306,x307,x308,x309,x310,x311
          ,x312,x313,x314,x315,x316,x317,x318,x319,x320,x321,x322,x323,x324
          ,x325,x326,x327,x328,x329,x330,x331,x332,x333,x334,x335,x336,x337
          ,x338,x339,x340,x341,x342,x343,x344,x345,x346,x347,x348,x349,x350
          ,x351,x352,x353,x354,x355,x356,x357,x358,x359,x360,x361,x362,x363
          ,x364,x365,x366,x367,x368,x369,x370,x371,x372,x373,x374,x375,x376
          ,x377,x378,x379,x380,x381,x382,x383,x384,x385,x386,x387,x388,x389
          ,x390,x391,x392,x393,x394,x395,x396,x397,x398,x399,x400,x401,x402
          ,x403,x404,x405,x406,x407,x408,x409,x410,x411,x412,x413,x414,x415
          ,x416,x417,x418,x419,x420,x421,x422,x423,x424,x425,x426,x427,x428
          ,x429,x430,x431,x432,x433,x434,x435,x436,x437,x438,x439,x440,x441
          ,x442,x443,x444,x445,x446,x447,x448,x449,x450,x451,x452,x453,x454
          ,x455,x456,x457,x458,x459,x460,x461,x462,x463,x464,x465,x466,x467
          ,x468,x469,x470,x471,x472,x473,x474,x475,x476,x477,x478,x479,x480
          ,x481,x482,x483,x484,x485,x486,x487,x488,x489,x490,x491,x492,x493
          ,x494,x495,x496,x497,x498,x499,x500,x501,x502,x503,x504,x505,objvar
          ,x507;

Positive Variables  x507;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19
          ,e20,e21,e22,e23,e24,e25,e26,e27,e28,e29,e30,e31,e32,e33,e34,e35,e36
          ,e37,e38,e39,e40,e41,e42,e43,e44,e45,e46,e47,e48,e49,e50,e51,e52,e53
          ,e54,e55,e56,e57,e58,e59,e60,e61,e62,e63,e64,e65,e66,e67,e68,e69,e70
          ,e71,e72,e73,e74,e75,e76,e77,e78,e79,e80,e81,e82,e83,e84,e85,e86,e87
          ,e88,e89,e90,e91,e92,e93,e94,e95,e96,e97,e98,e99,e100,e101,e102,e103
          ,e104,e105,e106,e107,e108,e109,e110,e111,e112,e113,e114,e115,e116
          ,e117,e118,e119,e120,e121,e122,e123,e124,e125,e126,e127,e128,e129
          ,e130,e131,e132,e133,e134,e135,e136,e137,e138,e139,e140,e141,e142
          ,e143,e144,e145,e146,e147,e148,e149,e150,e151,e152,e153,e154,e155
          ,e156,e157,e158,e159,e160,e161,e162,e163,e164,e165,e166,e167,e168
          ,e169,e170,e171,e172,e173,e174,e175,e176,e177,e178,e179,e180,e181
          ,e182,e183,e184,e185,e186,e187,e188,e189,e190,e191,e192,e193,e194
          ,e195,e196,e197,e198,e199,e200,e201,e202,e203,e204,e205,e206,e207
          ,e208,e209,e210,e211,e212,e213,e214,e215,e216,e217,e218,e219,e220
          ,e221,e222,e223,e224,e225,e226,e227,e228,e229,e230,e231,e232,e233
          ,e234,e235,e236,e237,e238,e239,e240,e241,e242,e243,e244,e245,e246
          ,e247,e248,e249,e250,e251,e252,e253,e254,e255,e256,e257,e258,e259
          ,e260,e261,e262,e263,e264,e265,e266,e267,e268,e269,e270,e271,e272
          ,e273,e274,e275,e276,e277,e278,e279,e280,e281,e282,e283,e284,e285
          ,e286,e287,e288,e289,e290,e291,e292,e293,e294,e295,e296,e297,e298
          ,e299,e300,e301,e302,e303,e304,e305,e306,e307,e308,e309,e310,e311
          ,e312,e313,e314,e315,e316,e317,e318,e319,e320,e321,e322,e323,e324
          ,e325,e326,e327,e328,e329,e330,e331,e332,e333,e334,e335,e336,e337
          ,e338,e339,e340,e341,e342,e343,e344,e345,e346,e347,e348,e349,e350
          ,e351,e352,e353,e354,e355,e356,e357,e358,e359,e360,e361,e362,e363
          ,e364,e365,e366,e367,e368,e369,e370,e371,e372,e373,e374,e375,e376
          ,e377,e378,e379,e380,e381,e382,e383,e384,e385,e386,e387,e388,e389
          ,e390,e391,e392,e393,e394,e395,e396,e397,e398,e399,e400,e401;


e1..    objvar - 100*x507 =E= 0;

e2.. -0.5*x507*(x304 + x305) - x102 + x103 =E= 0;

e3.. -0.5*x507*(x305 + x306) - x103 + x104 =E= 0;

e4.. -0.5*x507*(x306 + x307) - x104 + x105 =E= 0;

e5.. -0.5*x507*(x307 + x308) - x105 + x106 =E= 0;

e6.. -0.5*x507*(x308 + x309) - x106 + x107 =E= 0;

e7.. -0.5*x507*(x309 + x310) - x107 + x108 =E= 0;

e8.. -0.5*x507*(x310 + x311) - x108 + x109 =E= 0;

e9.. -0.5*x507*(x311 + x312) - x109 + x110 =E= 0;

e10.. -0.5*x507*(x312 + x313) - x110 + x111 =E= 0;

e11.. -0.5*x507*(x313 + x314) - x111 + x112 =E= 0;

e12.. -0.5*x507*(x314 + x315) - x112 + x113 =E= 0;

e13.. -0.5*x507*(x315 + x316) - x113 + x114 =E= 0;

e14.. -0.5*x507*(x316 + x317) - x114 + x115 =E= 0;

e15.. -0.5*x507*(x317 + x318) - x115 + x116 =E= 0;

e16.. -0.5*x507*(x318 + x319) - x116 + x117 =E= 0;

e17.. -0.5*x507*(x319 + x320) - x117 + x118 =E= 0;

e18.. -0.5*x507*(x320 + x321) - x118 + x119 =E= 0;

e19.. -0.5*x507*(x321 + x322) - x119 + x120 =E= 0;

e20.. -0.5*x507*(x322 + x323) - x120 + x121 =E= 0;

e21.. -0.5*x507*(x323 + x324) - x121 + x122 =E= 0;

e22.. -0.5*x507*(x324 + x325) - x122 + x123 =E= 0;

e23.. -0.5*x507*(x325 + x326) - x123 + x124 =E= 0;

e24.. -0.5*x507*(x326 + x327) - x124 + x125 =E= 0;

e25.. -0.5*x507*(x327 + x328) - x125 + x126 =E= 0;

e26.. -0.5*x507*(x328 + x329) - x126 + x127 =E= 0;

e27.. -0.5*x507*(x329 + x330) - x127 + x128 =E= 0;

e28.. -0.5*x507*(x330 + x331) - x128 + x129 =E= 0;

e29.. -0.5*x507*(x331 + x332) - x129 + x130 =E= 0;

e30.. -0.5*x507*(x332 + x333) - x130 + x131 =E= 0;

e31.. -0.5*x507*(x333 + x334) - x131 + x132 =E= 0;

e32.. -0.5*x507*(x334 + x335) - x132 + x133 =E= 0;

e33.. -0.5*x507*(x335 + x336) - x133 + x134 =E= 0;

e34.. -0.5*x507*(x336 + x337) - x134 + x135 =E= 0;

e35.. -0.5*x507*(x337 + x338) - x135 + x136 =E= 0;

e36.. -0.5*x507*(x338 + x339) - x136 + x137 =E= 0;

e37.. -0.5*x507*(x339 + x340) - x137 + x138 =E= 0;

e38.. -0.5*x507*(x340 + x341) - x138 + x139 =E= 0;

e39.. -0.5*x507*(x341 + x342) - x139 + x140 =E= 0;

e40.. -0.5*x507*(x342 + x343) - x140 + x141 =E= 0;

e41.. -0.5*x507*(x343 + x344) - x141 + x142 =E= 0;

e42.. -0.5*x507*(x344 + x345) - x142 + x143 =E= 0;

e43.. -0.5*x507*(x345 + x346) - x143 + x144 =E= 0;

e44.. -0.5*x507*(x346 + x347) - x144 + x145 =E= 0;

e45.. -0.5*x507*(x347 + x348) - x145 + x146 =E= 0;

e46.. -0.5*x507*(x348 + x349) - x146 + x147 =E= 0;

e47.. -0.5*x507*(x349 + x350) - x147 + x148 =E= 0;

e48.. -0.5*x507*(x350 + x351) - x148 + x149 =E= 0;

e49.. -0.5*x507*(x351 + x352) - x149 + x150 =E= 0;

e50.. -0.5*x507*(x352 + x353) - x150 + x151 =E= 0;

e51.. -0.5*x507*(x353 + x354) - x151 + x152 =E= 0;

e52.. -0.5*x507*(x354 + x355) - x152 + x153 =E= 0;

e53.. -0.5*x507*(x355 + x356) - x153 + x154 =E= 0;

e54.. -0.5*x507*(x356 + x357) - x154 + x155 =E= 0;

e55.. -0.5*x507*(x357 + x358) - x155 + x156 =E= 0;

e56.. -0.5*x507*(x358 + x359) - x156 + x157 =E= 0;

e57.. -0.5*x507*(x359 + x360) - x157 + x158 =E= 0;

e58.. -0.5*x507*(x360 + x361) - x158 + x159 =E= 0;

e59.. -0.5*x507*(x361 + x362) - x159 + x160 =E= 0;

e60.. -0.5*x507*(x362 + x363) - x160 + x161 =E= 0;

e61.. -0.5*x507*(x363 + x364) - x161 + x162 =E= 0;

e62.. -0.5*x507*(x364 + x365) - x162 + x163 =E= 0;

e63.. -0.5*x507*(x365 + x366) - x163 + x164 =E= 0;

e64.. -0.5*x507*(x366 + x367) - x164 + x165 =E= 0;

e65.. -0.5*x507*(x367 + x368) - x165 + x166 =E= 0;

e66.. -0.5*x507*(x368 + x369) - x166 + x167 =E= 0;

e67.. -0.5*x507*(x369 + x370) - x167 + x168 =E= 0;

e68.. -0.5*x507*(x370 + x371) - x168 + x169 =E= 0;

e69.. -0.5*x507*(x371 + x372) - x169 + x170 =E= 0;

e70.. -0.5*x507*(x372 + x373) - x170 + x171 =E= 0;

e71.. -0.5*x507*(x373 + x374) - x171 + x172 =E= 0;

e72.. -0.5*x507*(x374 + x375) - x172 + x173 =E= 0;

e73.. -0.5*x507*(x375 + x376) - x173 + x174 =E= 0;

e74.. -0.5*x507*(x376 + x377) - x174 + x175 =E= 0;

e75.. -0.5*x507*(x377 + x378) - x175 + x176 =E= 0;

e76.. -0.5*x507*(x378 + x379) - x176 + x177 =E= 0;

e77.. -0.5*x507*(x379 + x380) - x177 + x178 =E= 0;

e78.. -0.5*x507*(x380 + x381) - x178 + x179 =E= 0;

e79.. -0.5*x507*(x381 + x382) - x179 + x180 =E= 0;

e80.. -0.5*x507*(x382 + x383) - x180 + x181 =E= 0;

e81.. -0.5*x507*(x383 + x384) - x181 + x182 =E= 0;

e82.. -0.5*x507*(x384 + x385) - x182 + x183 =E= 0;

e83.. -0.5*x507*(x385 + x386) - x183 + x184 =E= 0;

e84.. -0.5*x507*(x386 + x387) - x184 + x185 =E= 0;

e85.. -0.5*x507*(x387 + x388) - x185 + x186 =E= 0;

e86.. -0.5*x507*(x388 + x389) - x186 + x187 =E= 0;

e87.. -0.5*x507*(x389 + x390) - x187 + x188 =E= 0;

e88.. -0.5*x507*(x390 + x391) - x188 + x189 =E= 0;

e89.. -0.5*x507*(x391 + x392) - x189 + x190 =E= 0;

e90.. -0.5*x507*(x392 + x393) - x190 + x191 =E= 0;

e91.. -0.5*x507*(x393 + x394) - x191 + x192 =E= 0;

e92.. -0.5*x507*(x394 + x395) - x192 + x193 =E= 0;

e93.. -0.5*x507*(x395 + x396) - x193 + x194 =E= 0;

e94.. -0.5*x507*(x396 + x397) - x194 + x195 =E= 0;

e95.. -0.5*x507*(x397 + x398) - x195 + x196 =E= 0;

e96.. -0.5*x507*(x398 + x399) - x196 + x197 =E= 0;

e97.. -0.5*x507*(x399 + x400) - x197 + x198 =E= 0;

e98.. -0.5*x507*(x400 + x401) - x198 + x199 =E= 0;

e99.. -0.5*x507*(x401 + x402) - x199 + x200 =E= 0;

e100.. -0.5*x507*(x402 + x403) - x200 + x201 =E= 0;

e101.. -0.5*x507*(x403 + x404) - x201 + x202 =E= 0;

e102.. -0.5*x507*(x405 + x406) - x203 + x204 =E= 0;

e103.. -0.5*x507*(x406 + x407) - x204 + x205 =E= 0;

e104.. -0.5*x507*(x407 + x408) - x205 + x206 =E= 0;

e105.. -0.5*x507*(x408 + x409) - x206 + x207 =E= 0;

e106.. -0.5*x507*(x409 + x410) - x207 + x208 =E= 0;

e107.. -0.5*x507*(x410 + x411) - x208 + x209 =E= 0;

e108.. -0.5*x507*(x411 + x412) - x209 + x210 =E= 0;

e109.. -0.5*x507*(x412 + x413) - x210 + x211 =E= 0;

e110.. -0.5*x507*(x413 + x414) - x211 + x212 =E= 0;

e111.. -0.5*x507*(x414 + x415) - x212 + x213 =E= 0;

e112.. -0.5*x507*(x415 + x416) - x213 + x214 =E= 0;

e113.. -0.5*x507*(x416 + x417) - x214 + x215 =E= 0;

e114.. -0.5*x507*(x417 + x418) - x215 + x216 =E= 0;

e115.. -0.5*x507*(x418 + x419) - x216 + x217 =E= 0;

e116.. -0.5*x507*(x419 + x420) - x217 + x218 =E= 0;

e117.. -0.5*x507*(x420 + x421) - x218 + x219 =E= 0;

e118.. -0.5*x507*(x421 + x422) - x219 + x220 =E= 0;

e119.. -0.5*x507*(x422 + x423) - x220 + x221 =E= 0;

e120.. -0.5*x507*(x423 + x424) - x221 + x222 =E= 0;

e121.. -0.5*x507*(x424 + x425) - x222 + x223 =E= 0;

e122.. -0.5*x507*(x425 + x426) - x223 + x224 =E= 0;

e123.. -0.5*x507*(x426 + x427) - x224 + x225 =E= 0;

e124.. -0.5*x507*(x427 + x428) - x225 + x226 =E= 0;

e125.. -0.5*x507*(x428 + x429) - x226 + x227 =E= 0;

e126.. -0.5*x507*(x429 + x430) - x227 + x228 =E= 0;

e127.. -0.5*x507*(x430 + x431) - x228 + x229 =E= 0;

e128.. -0.5*x507*(x431 + x432) - x229 + x230 =E= 0;

e129.. -0.5*x507*(x432 + x433) - x230 + x231 =E= 0;

e130.. -0.5*x507*(x433 + x434) - x231 + x232 =E= 0;

e131.. -0.5*x507*(x434 + x435) - x232 + x233 =E= 0;

e132.. -0.5*x507*(x435 + x436) - x233 + x234 =E= 0;

e133.. -0.5*x507*(x436 + x437) - x234 + x235 =E= 0;

e134.. -0.5*x507*(x437 + x438) - x235 + x236 =E= 0;

e135.. -0.5*x507*(x438 + x439) - x236 + x237 =E= 0;

e136.. -0.5*x507*(x439 + x440) - x237 + x238 =E= 0;

e137.. -0.5*x507*(x440 + x441) - x238 + x239 =E= 0;

e138.. -0.5*x507*(x441 + x442) - x239 + x240 =E= 0;

e139.. -0.5*x507*(x442 + x443) - x240 + x241 =E= 0;

e140.. -0.5*x507*(x443 + x444) - x241 + x242 =E= 0;

e141.. -0.5*x507*(x444 + x445) - x242 + x243 =E= 0;

e142.. -0.5*x507*(x445 + x446) - x243 + x244 =E= 0;

e143.. -0.5*x507*(x446 + x447) - x244 + x245 =E= 0;

e144.. -0.5*x507*(x447 + x448) - x245 + x246 =E= 0;

e145.. -0.5*x507*(x448 + x449) - x246 + x247 =E= 0;

e146.. -0.5*x507*(x449 + x450) - x247 + x248 =E= 0;

e147.. -0.5*x507*(x450 + x451) - x248 + x249 =E= 0;

e148.. -0.5*x507*(x451 + x452) - x249 + x250 =E= 0;

e149.. -0.5*x507*(x452 + x453) - x250 + x251 =E= 0;

e150.. -0.5*x507*(x453 + x454) - x251 + x252 =E= 0;

e151.. -0.5*x507*(x454 + x455) - x252 + x253 =E= 0;

e152.. -0.5*x507*(x455 + x456) - x253 + x254 =E= 0;

e153.. -0.5*x507*(x456 + x457) - x254 + x255 =E= 0;

e154.. -0.5*x507*(x457 + x458) - x255 + x256 =E= 0;

e155.. -0.5*x507*(x458 + x459) - x256 + x257 =E= 0;

e156.. -0.5*x507*(x459 + x460) - x257 + x258 =E= 0;

e157.. -0.5*x507*(x460 + x461) - x258 + x259 =E= 0;

e158.. -0.5*x507*(x461 + x462) - x259 + x260 =E= 0;

e159.. -0.5*x507*(x462 + x463) - x260 + x261 =E= 0;

e160.. -0.5*x507*(x463 + x464) - x261 + x262 =E= 0;

e161.. -0.5*x507*(x464 + x465) - x262 + x263 =E= 0;

e162.. -0.5*x507*(x465 + x466) - x263 + x264 =E= 0;

e163.. -0.5*x507*(x466 + x467) - x264 + x265 =E= 0;

e164.. -0.5*x507*(x467 + x468) - x265 + x266 =E= 0;

e165.. -0.5*x507*(x468 + x469) - x266 + x267 =E= 0;

e166.. -0.5*x507*(x469 + x470) - x267 + x268 =E= 0;

e167.. -0.5*x507*(x470 + x471) - x268 + x269 =E= 0;

e168.. -0.5*x507*(x471 + x472) - x269 + x270 =E= 0;

e169.. -0.5*x507*(x472 + x473) - x270 + x271 =E= 0;

e170.. -0.5*x507*(x473 + x474) - x271 + x272 =E= 0;

e171.. -0.5*x507*(x474 + x475) - x272 + x273 =E= 0;

e172.. -0.5*x507*(x475 + x476) - x273 + x274 =E= 0;

e173.. -0.5*x507*(x476 + x477) - x274 + x275 =E= 0;

e174.. -0.5*x507*(x477 + x478) - x275 + x276 =E= 0;

e175.. -0.5*x507*(x478 + x479) - x276 + x277 =E= 0;

e176.. -0.5*x507*(x479 + x480) - x277 + x278 =E= 0;

e177.. -0.5*x507*(x480 + x481) - x278 + x279 =E= 0;

e178.. -0.5*x507*(x481 + x482) - x279 + x280 =E= 0;

e179.. -0.5*x507*(x482 + x483) - x280 + x281 =E= 0;

e180.. -0.5*x507*(x483 + x484) - x281 + x282 =E= 0;

e181.. -0.5*x507*(x484 + x485) - x282 + x283 =E= 0;

e182.. -0.5*x507*(x485 + x486) - x283 + x284 =E= 0;

e183.. -0.5*x507*(x486 + x487) - x284 + x285 =E= 0;

e184.. -0.5*x507*(x487 + x488) - x285 + x286 =E= 0;

e185.. -0.5*x507*(x488 + x489) - x286 + x287 =E= 0;

e186.. -0.5*x507*(x489 + x490) - x287 + x288 =E= 0;

e187.. -0.5*x507*(x490 + x491) - x288 + x289 =E= 0;

e188.. -0.5*x507*(x491 + x492) - x289 + x290 =E= 0;

e189.. -0.5*x507*(x492 + x493) - x290 + x291 =E= 0;

e190.. -0.5*x507*(x493 + x494) - x291 + x292 =E= 0;

e191.. -0.5*x507*(x494 + x495) - x292 + x293 =E= 0;

e192.. -0.5*x507*(x495 + x496) - x293 + x294 =E= 0;

e193.. -0.5*x507*(x496 + x497) - x294 + x295 =E= 0;

e194.. -0.5*x507*(x497 + x498) - x295 + x296 =E= 0;

e195.. -0.5*x507*(x498 + x499) - x296 + x297 =E= 0;

e196.. -0.5*x507*(x499 + x500) - x297 + x298 =E= 0;

e197.. -0.5*x507*(x500 + x501) - x298 + x299 =E= 0;

e198.. -0.5*x507*(x501 + x502) - x299 + x300 =E= 0;

e199.. -0.5*x507*(x502 + x503) - x300 + x301 =E= 0;

e200.. -0.5*x507*(x503 + x504) - x301 + x302 =E= 0;

e201.. -0.5*x507*(x504 + x505) - x302 + x303 =E= 0;

e202.. -0.5*(100*cos(x1) + 100*cos(x2))*x507 - x304 + x305 =E= 0;

e203.. -0.5*(100*cos(x2) + 100*cos(x3))*x507 - x305 + x306 =E= 0;

e204.. -0.5*(100*cos(x3) + 100*cos(x4))*x507 - x306 + x307 =E= 0;

e205.. -0.5*(100*cos(x4) + 100*cos(x5))*x507 - x307 + x308 =E= 0;

e206.. -0.5*(100*cos(x5) + 100*cos(x6))*x507 - x308 + x309 =E= 0;

e207.. -0.5*(100*cos(x6) + 100*cos(x7))*x507 - x309 + x310 =E= 0;

e208.. -0.5*(100*cos(x7) + 100*cos(x8))*x507 - x310 + x311 =E= 0;

e209.. -0.5*(100*cos(x8) + 100*cos(x9))*x507 - x311 + x312 =E= 0;

e210.. -0.5*(100*cos(x9) + 100*cos(x10))*x507 - x312 + x313 =E= 0;

e211.. -0.5*(100*cos(x10) + 100*cos(x11))*x507 - x313 + x314 =E= 0;

e212.. -0.5*(100*cos(x11) + 100*cos(x12))*x507 - x314 + x315 =E= 0;

e213.. -0.5*(100*cos(x12) + 100*cos(x13))*x507 - x315 + x316 =E= 0;

e214.. -0.5*(100*cos(x13) + 100*cos(x14))*x507 - x316 + x317 =E= 0;

e215.. -0.5*(100*cos(x14) + 100*cos(x15))*x507 - x317 + x318 =E= 0;

e216.. -0.5*(100*cos(x15) + 100*cos(x16))*x507 - x318 + x319 =E= 0;

e217.. -0.5*(100*cos(x16) + 100*cos(x17))*x507 - x319 + x320 =E= 0;

e218.. -0.5*(100*cos(x17) + 100*cos(x18))*x507 - x320 + x321 =E= 0;

e219.. -0.5*(100*cos(x18) + 100*cos(x19))*x507 - x321 + x322 =E= 0;

e220.. -0.5*(100*cos(x19) + 100*cos(x20))*x507 - x322 + x323 =E= 0;

e221.. -0.5*(100*cos(x20) + 100*cos(x21))*x507 - x323 + x324 =E= 0;

e222.. -0.5*(100*cos(x21) + 100*cos(x22))*x507 - x324 + x325 =E= 0;

e223.. -0.5*(100*cos(x22) + 100*cos(x23))*x507 - x325 + x326 =E= 0;

e224.. -0.5*(100*cos(x23) + 100*cos(x24))*x507 - x326 + x327 =E= 0;

e225.. -0.5*(100*cos(x24) + 100*cos(x25))*x507 - x327 + x328 =E= 0;

e226.. -0.5*(100*cos(x25) + 100*cos(x26))*x507 - x328 + x329 =E= 0;

e227.. -0.5*(100*cos(x26) + 100*cos(x27))*x507 - x329 + x330 =E= 0;

e228.. -0.5*(100*cos(x27) + 100*cos(x28))*x507 - x330 + x331 =E= 0;

e229.. -0.5*(100*cos(x28) + 100*cos(x29))*x507 - x331 + x332 =E= 0;

e230.. -0.5*(100*cos(x29) + 100*cos(x30))*x507 - x332 + x333 =E= 0;

e231.. -0.5*(100*cos(x30) + 100*cos(x31))*x507 - x333 + x334 =E= 0;

e232.. -0.5*(100*cos(x31) + 100*cos(x32))*x507 - x334 + x335 =E= 0;

e233.. -0.5*(100*cos(x32) + 100*cos(x33))*x507 - x335 + x336 =E= 0;

e234.. -0.5*(100*cos(x33) + 100*cos(x34))*x507 - x336 + x337 =E= 0;

e235.. -0.5*(100*cos(x34) + 100*cos(x35))*x507 - x337 + x338 =E= 0;

e236.. -0.5*(100*cos(x35) + 100*cos(x36))*x507 - x338 + x339 =E= 0;

e237.. -0.5*(100*cos(x36) + 100*cos(x37))*x507 - x339 + x340 =E= 0;

e238.. -0.5*(100*cos(x37) + 100*cos(x38))*x507 - x340 + x341 =E= 0;

e239.. -0.5*(100*cos(x38) + 100*cos(x39))*x507 - x341 + x342 =E= 0;

e240.. -0.5*(100*cos(x39) + 100*cos(x40))*x507 - x342 + x343 =E= 0;

e241.. -0.5*(100*cos(x40) + 100*cos(x41))*x507 - x343 + x344 =E= 0;

e242.. -0.5*(100*cos(x41) + 100*cos(x42))*x507 - x344 + x345 =E= 0;

e243.. -0.5*(100*cos(x42) + 100*cos(x43))*x507 - x345 + x346 =E= 0;

e244.. -0.5*(100*cos(x43) + 100*cos(x44))*x507 - x346 + x347 =E= 0;

e245.. -0.5*(100*cos(x44) + 100*cos(x45))*x507 - x347 + x348 =E= 0;

e246.. -0.5*(100*cos(x45) + 100*cos(x46))*x507 - x348 + x349 =E= 0;

e247.. -0.5*(100*cos(x46) + 100*cos(x47))*x507 - x349 + x350 =E= 0;

e248.. -0.5*(100*cos(x47) + 100*cos(x48))*x507 - x350 + x351 =E= 0;

e249.. -0.5*(100*cos(x48) + 100*cos(x49))*x507 - x351 + x352 =E= 0;

e250.. -0.5*(100*cos(x49) + 100*cos(x50))*x507 - x352 + x353 =E= 0;

e251.. -0.5*(100*cos(x50) + 100*cos(x51))*x507 - x353 + x354 =E= 0;

e252.. -0.5*(100*cos(x51) + 100*cos(x52))*x507 - x354 + x355 =E= 0;

e253.. -0.5*(100*cos(x52) + 100*cos(x53))*x507 - x355 + x356 =E= 0;

e254.. -0.5*(100*cos(x53) + 100*cos(x54))*x507 - x356 + x357 =E= 0;

e255.. -0.5*(100*cos(x54) + 100*cos(x55))*x507 - x357 + x358 =E= 0;

e256.. -0.5*(100*cos(x55) + 100*cos(x56))*x507 - x358 + x359 =E= 0;

e257.. -0.5*(100*cos(x56) + 100*cos(x57))*x507 - x359 + x360 =E= 0;

e258.. -0.5*(100*cos(x57) + 100*cos(x58))*x507 - x360 + x361 =E= 0;

e259.. -0.5*(100*cos(x58) + 100*cos(x59))*x507 - x361 + x362 =E= 0;

e260.. -0.5*(100*cos(x59) + 100*cos(x60))*x507 - x362 + x363 =E= 0;

e261.. -0.5*(100*cos(x60) + 100*cos(x61))*x507 - x363 + x364 =E= 0;

e262.. -0.5*(100*cos(x61) + 100*cos(x62))*x507 - x364 + x365 =E= 0;

e263.. -0.5*(100*cos(x62) + 100*cos(x63))*x507 - x365 + x366 =E= 0;

e264.. -0.5*(100*cos(x63) + 100*cos(x64))*x507 - x366 + x367 =E= 0;

e265.. -0.5*(100*cos(x64) + 100*cos(x65))*x507 - x367 + x368 =E= 0;

e266.. -0.5*(100*cos(x65) + 100*cos(x66))*x507 - x368 + x369 =E= 0;

e267.. -0.5*(100*cos(x66) + 100*cos(x67))*x507 - x369 + x370 =E= 0;

e268.. -0.5*(100*cos(x67) + 100*cos(x68))*x507 - x370 + x371 =E= 0;

e269.. -0.5*(100*cos(x68) + 100*cos(x69))*x507 - x371 + x372 =E= 0;

e270.. -0.5*(100*cos(x69) + 100*cos(x70))*x507 - x372 + x373 =E= 0;

e271.. -0.5*(100*cos(x70) + 100*cos(x71))*x507 - x373 + x374 =E= 0;

e272.. -0.5*(100*cos(x71) + 100*cos(x72))*x507 - x374 + x375 =E= 0;

e273.. -0.5*(100*cos(x72) + 100*cos(x73))*x507 - x375 + x376 =E= 0;

e274.. -0.5*(100*cos(x73) + 100*cos(x74))*x507 - x376 + x377 =E= 0;

e275.. -0.5*(100*cos(x74) + 100*cos(x75))*x507 - x377 + x378 =E= 0;

e276.. -0.5*(100*cos(x75) + 100*cos(x76))*x507 - x378 + x379 =E= 0;

e277.. -0.5*(100*cos(x76) + 100*cos(x77))*x507 - x379 + x380 =E= 0;

e278.. -0.5*(100*cos(x77) + 100*cos(x78))*x507 - x380 + x381 =E= 0;

e279.. -0.5*(100*cos(x78) + 100*cos(x79))*x507 - x381 + x382 =E= 0;

e280.. -0.5*(100*cos(x79) + 100*cos(x80))*x507 - x382 + x383 =E= 0;

e281.. -0.5*(100*cos(x80) + 100*cos(x81))*x507 - x383 + x384 =E= 0;

e282.. -0.5*(100*cos(x81) + 100*cos(x82))*x507 - x384 + x385 =E= 0;

e283.. -0.5*(100*cos(x82) + 100*cos(x83))*x507 - x385 + x386 =E= 0;

e284.. -0.5*(100*cos(x83) + 100*cos(x84))*x507 - x386 + x387 =E= 0;

e285.. -0.5*(100*cos(x84) + 100*cos(x85))*x507 - x387 + x388 =E= 0;

e286.. -0.5*(100*cos(x85) + 100*cos(x86))*x507 - x388 + x389 =E= 0;

e287.. -0.5*(100*cos(x86) + 100*cos(x87))*x507 - x389 + x390 =E= 0;

e288.. -0.5*(100*cos(x87) + 100*cos(x88))*x507 - x390 + x391 =E= 0;

e289.. -0.5*(100*cos(x88) + 100*cos(x89))*x507 - x391 + x392 =E= 0;

e290.. -0.5*(100*cos(x89) + 100*cos(x90))*x507 - x392 + x393 =E= 0;

e291.. -0.5*(100*cos(x90) + 100*cos(x91))*x507 - x393 + x394 =E= 0;

e292.. -0.5*(100*cos(x91) + 100*cos(x92))*x507 - x394 + x395 =E= 0;

e293.. -0.5*(100*cos(x92) + 100*cos(x93))*x507 - x395 + x396 =E= 0;

e294.. -0.5*(100*cos(x93) + 100*cos(x94))*x507 - x396 + x397 =E= 0;

e295.. -0.5*(100*cos(x94) + 100*cos(x95))*x507 - x397 + x398 =E= 0;

e296.. -0.5*(100*cos(x95) + 100*cos(x96))*x507 - x398 + x399 =E= 0;

e297.. -0.5*(100*cos(x96) + 100*cos(x97))*x507 - x399 + x400 =E= 0;

e298.. -0.5*(100*cos(x97) + 100*cos(x98))*x507 - x400 + x401 =E= 0;

e299.. -0.5*(100*cos(x98) + 100*cos(x99))*x507 - x401 + x402 =E= 0;

e300.. -0.5*(100*cos(x99) + 100*cos(x100))*x507 - x402 + x403 =E= 0;

e301.. -0.5*(100*cos(x100) + 100*cos(x101))*x507 - x403 + x404 =E= 0;

e302.. -0.5*(100*sin(x1) + 100*sin(x2))*x507 - x405 + x406 =E= 0;

e303.. -0.5*(100*sin(x2) + 100*sin(x3))*x507 - x406 + x407 =E= 0;

e304.. -0.5*(100*sin(x3) + 100*sin(x4))*x507 - x407 + x408 =E= 0;

e305.. -0.5*(100*sin(x4) + 100*sin(x5))*x507 - x408 + x409 =E= 0;

e306.. -0.5*(100*sin(x5) + 100*sin(x6))*x507 - x409 + x410 =E= 0;

e307.. -0.5*(100*sin(x6) + 100*sin(x7))*x507 - x410 + x411 =E= 0;

e308.. -0.5*(100*sin(x7) + 100*sin(x8))*x507 - x411 + x412 =E= 0;

e309.. -0.5*(100*sin(x8) + 100*sin(x9))*x507 - x412 + x413 =E= 0;

e310.. -0.5*(100*sin(x9) + 100*sin(x10))*x507 - x413 + x414 =E= 0;

e311.. -0.5*(100*sin(x10) + 100*sin(x11))*x507 - x414 + x415 =E= 0;

e312.. -0.5*(100*sin(x11) + 100*sin(x12))*x507 - x415 + x416 =E= 0;

e313.. -0.5*(100*sin(x12) + 100*sin(x13))*x507 - x416 + x417 =E= 0;

e314.. -0.5*(100*sin(x13) + 100*sin(x14))*x507 - x417 + x418 =E= 0;

e315.. -0.5*(100*sin(x14) + 100*sin(x15))*x507 - x418 + x419 =E= 0;

e316.. -0.5*(100*sin(x15) + 100*sin(x16))*x507 - x419 + x420 =E= 0;

e317.. -0.5*(100*sin(x16) + 100*sin(x17))*x507 - x420 + x421 =E= 0;

e318.. -0.5*(100*sin(x17) + 100*sin(x18))*x507 - x421 + x422 =E= 0;

e319.. -0.5*(100*sin(x18) + 100*sin(x19))*x507 - x422 + x423 =E= 0;

e320.. -0.5*(100*sin(x19) + 100*sin(x20))*x507 - x423 + x424 =E= 0;

e321.. -0.5*(100*sin(x20) + 100*sin(x21))*x507 - x424 + x425 =E= 0;

e322.. -0.5*(100*sin(x21) + 100*sin(x22))*x507 - x425 + x426 =E= 0;

e323.. -0.5*(100*sin(x22) + 100*sin(x23))*x507 - x426 + x427 =E= 0;

e324.. -0.5*(100*sin(x23) + 100*sin(x24))*x507 - x427 + x428 =E= 0;

e325.. -0.5*(100*sin(x24) + 100*sin(x25))*x507 - x428 + x429 =E= 0;

e326.. -0.5*(100*sin(x25) + 100*sin(x26))*x507 - x429 + x430 =E= 0;

e327.. -0.5*(100*sin(x26) + 100*sin(x27))*x507 - x430 + x431 =E= 0;

e328.. -0.5*(100*sin(x27) + 100*sin(x28))*x507 - x431 + x432 =E= 0;

e329.. -0.5*(100*sin(x28) + 100*sin(x29))*x507 - x432 + x433 =E= 0;

e330.. -0.5*(100*sin(x29) + 100*sin(x30))*x507 - x433 + x434 =E= 0;

e331.. -0.5*(100*sin(x30) + 100*sin(x31))*x507 - x434 + x435 =E= 0;

e332.. -0.5*(100*sin(x31) + 100*sin(x32))*x507 - x435 + x436 =E= 0;

e333.. -0.5*(100*sin(x32) + 100*sin(x33))*x507 - x436 + x437 =E= 0;

e334.. -0.5*(100*sin(x33) + 100*sin(x34))*x507 - x437 + x438 =E= 0;

e335.. -0.5*(100*sin(x34) + 100*sin(x35))*x507 - x438 + x439 =E= 0;

e336.. -0.5*(100*sin(x35) + 100*sin(x36))*x507 - x439 + x440 =E= 0;

e337.. -0.5*(100*sin(x36) + 100*sin(x37))*x507 - x440 + x441 =E= 0;

e338.. -0.5*(100*sin(x37) + 100*sin(x38))*x507 - x441 + x442 =E= 0;

e339.. -0.5*(100*sin(x38) + 100*sin(x39))*x507 - x442 + x443 =E= 0;

e340.. -0.5*(100*sin(x39) + 100*sin(x40))*x507 - x443 + x444 =E= 0;

e341.. -0.5*(100*sin(x40) + 100*sin(x41))*x507 - x444 + x445 =E= 0;

e342.. -0.5*(100*sin(x41) + 100*sin(x42))*x507 - x445 + x446 =E= 0;

e343.. -0.5*(100*sin(x42) + 100*sin(x43))*x507 - x446 + x447 =E= 0;

e344.. -0.5*(100*sin(x43) + 100*sin(x44))*x507 - x447 + x448 =E= 0;

e345.. -0.5*(100*sin(x44) + 100*sin(x45))*x507 - x448 + x449 =E= 0;

e346.. -0.5*(100*sin(x45) + 100*sin(x46))*x507 - x449 + x450 =E= 0;

e347.. -0.5*(100*sin(x46) + 100*sin(x47))*x507 - x450 + x451 =E= 0;

e348.. -0.5*(100*sin(x47) + 100*sin(x48))*x507 - x451 + x452 =E= 0;

e349.. -0.5*(100*sin(x48) + 100*sin(x49))*x507 - x452 + x453 =E= 0;

e350.. -0.5*(100*sin(x49) + 100*sin(x50))*x507 - x453 + x454 =E= 0;

e351.. -0.5*(100*sin(x50) + 100*sin(x51))*x507 - x454 + x455 =E= 0;

e352.. -0.5*(100*sin(x51) + 100*sin(x52))*x507 - x455 + x456 =E= 0;

e353.. -0.5*(100*sin(x52) + 100*sin(x53))*x507 - x456 + x457 =E= 0;

e354.. -0.5*(100*sin(x53) + 100*sin(x54))*x507 - x457 + x458 =E= 0;

e355.. -0.5*(100*sin(x54) + 100*sin(x55))*x507 - x458 + x459 =E= 0;

e356.. -0.5*(100*sin(x55) + 100*sin(x56))*x507 - x459 + x460 =E= 0;

e357.. -0.5*(100*sin(x56) + 100*sin(x57))*x507 - x460 + x461 =E= 0;

e358.. -0.5*(100*sin(x57) + 100*sin(x58))*x507 - x461 + x462 =E= 0;

e359.. -0.5*(100*sin(x58) + 100*sin(x59))*x507 - x462 + x463 =E= 0;

e360.. -0.5*(100*sin(x59) + 100*sin(x60))*x507 - x463 + x464 =E= 0;

e361.. -0.5*(100*sin(x60) + 100*sin(x61))*x507 - x464 + x465 =E= 0;

e362.. -0.5*(100*sin(x61) + 100*sin(x62))*x507 - x465 + x466 =E= 0;

e363.. -0.5*(100*sin(x62) + 100*sin(x63))*x507 - x466 + x467 =E= 0;

e364.. -0.5*(100*sin(x63) + 100*sin(x64))*x507 - x467 + x468 =E= 0;

e365.. -0.5*(100*sin(x64) + 100*sin(x65))*x507 - x468 + x469 =E= 0;

e366.. -0.5*(100*sin(x65) + 100*sin(x66))*x507 - x469 + x470 =E= 0;

e367.. -0.5*(100*sin(x66) + 100*sin(x67))*x507 - x470 + x471 =E= 0;

e368.. -0.5*(100*sin(x67) + 100*sin(x68))*x507 - x471 + x472 =E= 0;

e369.. -0.5*(100*sin(x68) + 100*sin(x69))*x507 - x472 + x473 =E= 0;

e370.. -0.5*(100*sin(x69) + 100*sin(x70))*x507 - x473 + x474 =E= 0;

e371.. -0.5*(100*sin(x70) + 100*sin(x71))*x507 - x474 + x475 =E= 0;

e372.. -0.5*(100*sin(x71) + 100*sin(x72))*x507 - x475 + x476 =E= 0;

e373.. -0.5*(100*sin(x72) + 100*sin(x73))*x507 - x476 + x477 =E= 0;

e374.. -0.5*(100*sin(x73) + 100*sin(x74))*x507 - x477 + x478 =E= 0;

e375.. -0.5*(100*sin(x74) + 100*sin(x75))*x507 - x478 + x479 =E= 0;

e376.. -0.5*(100*sin(x75) + 100*sin(x76))*x507 - x479 + x480 =E= 0;

e377.. -0.5*(100*sin(x76) + 100*sin(x77))*x507 - x480 + x481 =E= 0;

e378.. -0.5*(100*sin(x77) + 100*sin(x78))*x507 - x481 + x482 =E= 0;

e379.. -0.5*(100*sin(x78) + 100*sin(x79))*x507 - x482 + x483 =E= 0;

e380.. -0.5*(100*sin(x79) + 100*sin(x80))*x507 - x483 + x484 =E= 0;

e381.. -0.5*(100*sin(x80) + 100*sin(x81))*x507 - x484 + x485 =E= 0;

e382.. -0.5*(100*sin(x81) + 100*sin(x82))*x507 - x485 + x486 =E= 0;

e383.. -0.5*(100*sin(x82) + 100*sin(x83))*x507 - x486 + x487 =E= 0;

e384.. -0.5*(100*sin(x83) + 100*sin(x84))*x507 - x487 + x488 =E= 0;

e385.. -0.5*(100*sin(x84) + 100*sin(x85))*x507 - x488 + x489 =E= 0;

e386.. -0.5*(100*sin(x85) + 100*sin(x86))*x507 - x489 + x490 =E= 0;

e387.. -0.5*(100*sin(x86) + 100*sin(x87))*x507 - x490 + x491 =E= 0;

e388.. -0.5*(100*sin(x87) + 100*sin(x88))*x507 - x491 + x492 =E= 0;

e389.. -0.5*(100*sin(x88) + 100*sin(x89))*x507 - x492 + x493 =E= 0;

e390.. -0.5*(100*sin(x89) + 100*sin(x90))*x507 - x493 + x494 =E= 0;

e391.. -0.5*(100*sin(x90) + 100*sin(x91))*x507 - x494 + x495 =E= 0;

e392.. -0.5*(100*sin(x91) + 100*sin(x92))*x507 - x495 + x496 =E= 0;

e393.. -0.5*(100*sin(x92) + 100*sin(x93))*x507 - x496 + x497 =E= 0;

e394.. -0.5*(100*sin(x93) + 100*sin(x94))*x507 - x497 + x498 =E= 0;

e395.. -0.5*(100*sin(x94) + 100*sin(x95))*x507 - x498 + x499 =E= 0;

e396.. -0.5*(100*sin(x95) + 100*sin(x96))*x507 - x499 + x500 =E= 0;

e397.. -0.5*(100*sin(x96) + 100*sin(x97))*x507 - x500 + x501 =E= 0;

e398.. -0.5*(100*sin(x97) + 100*sin(x98))*x507 - x501 + x502 =E= 0;

e399.. -0.5*(100*sin(x98) + 100*sin(x99))*x507 - x502 + x503 =E= 0;

e400.. -0.5*(100*sin(x99) + 100*sin(x100))*x507 - x503 + x504 =E= 0;

e401.. -0.5*(100*sin(x100) + 100*sin(x101))*x507 - x504 + x505 =E= 0;

* set non-default bounds
x1.lo = -1.5707963267949; x1.up = 1.5707963267949;
x2.lo = -1.5707963267949; x2.up = 1.5707963267949;
x3.lo = -1.5707963267949; x3.up = 1.5707963267949;
x4.lo = -1.5707963267949; x4.up = 1.5707963267949;
x5.lo = -1.5707963267949; x5.up = 1.5707963267949;
x6.lo = -1.5707963267949; x6.up = 1.5707963267949;
x7.lo = -1.5707963267949; x7.up = 1.5707963267949;
x8.lo = -1.5707963267949; x8.up = 1.5707963267949;
x9.lo = -1.5707963267949; x9.up = 1.5707963267949;
x10.lo = -1.5707963267949; x10.up = 1.5707963267949;
x11.lo = -1.5707963267949; x11.up = 1.5707963267949;
x12.lo = -1.5707963267949; x12.up = 1.5707963267949;
x13.lo = -1.5707963267949; x13.up = 1.5707963267949;
x14.lo = -1.5707963267949; x14.up = 1.5707963267949;
x15.lo = -1.5707963267949; x15.up = 1.5707963267949;
x16.lo = -1.5707963267949; x16.up = 1.5707963267949;
x17.lo = -1.5707963267949; x17.up = 1.5707963267949;
x18.lo = -1.5707963267949; x18.up = 1.5707963267949;
x19.lo = -1.5707963267949; x19.up = 1.5707963267949;
x20.lo = -1.5707963267949; x20.up = 1.5707963267949;
x21.lo = -1.5707963267949; x21.up = 1.5707963267949;
x22.lo = -1.5707963267949; x22.up = 1.5707963267949;
x23.lo = -1.5707963267949; x23.up = 1.5707963267949;
x24.lo = -1.5707963267949; x24.up = 1.5707963267949;
x25.lo = -1.5707963267949; x25.up = 1.5707963267949;
x26.lo = -1.5707963267949; x26.up = 1.5707963267949;
x27.lo = -1.5707963267949; x27.up = 1.5707963267949;
x28.lo = -1.5707963267949; x28.up = 1.5707963267949;
x29.lo = -1.5707963267949; x29.up = 1.5707963267949;
x30.lo = -1.5707963267949; x30.up = 1.5707963267949;
x31.lo = -1.5707963267949; x31.up = 1.5707963267949;
x32.lo = -1.5707963267949; x32.up = 1.5707963267949;
x33.lo = -1.5707963267949; x33.up = 1.5707963267949;
x34.lo = -1.5707963267949; x34.up = 1.5707963267949;
x35.lo = -1.5707963267949; x35.up = 1.5707963267949;
x36.lo = -1.5707963267949; x36.up = 1.5707963267949;
x37.lo = -1.5707963267949; x37.up = 1.5707963267949;
x38.lo = -1.5707963267949; x38.up = 1.5707963267949;
x39.lo = -1.5707963267949; x39.up = 1.5707963267949;
x40.lo = -1.5707963267949; x40.up = 1.5707963267949;
x41.lo = -1.5707963267949; x41.up = 1.5707963267949;
x42.lo = -1.5707963267949; x42.up = 1.5707963267949;
x43.lo = -1.5707963267949; x43.up = 1.5707963267949;
x44.lo = -1.5707963267949; x44.up = 1.5707963267949;
x45.lo = -1.5707963267949; x45.up = 1.5707963267949;
x46.lo = -1.5707963267949; x46.up = 1.5707963267949;
x47.lo = -1.5707963267949; x47.up = 1.5707963267949;
x48.lo = -1.5707963267949; x48.up = 1.5707963267949;
x49.lo = -1.5707963267949; x49.up = 1.5707963267949;
x50.lo = -1.5707963267949; x50.up = 1.5707963267949;
x51.lo = -1.5707963267949; x51.up = 1.5707963267949;
x52.lo = -1.5707963267949; x52.up = 1.5707963267949;
x53.lo = -1.5707963267949; x53.up = 1.5707963267949;
x54.lo = -1.5707963267949; x54.up = 1.5707963267949;
x55.lo = -1.5707963267949; x55.up = 1.5707963267949;
x56.lo = -1.5707963267949; x56.up = 1.5707963267949;
x57.lo = -1.5707963267949; x57.up = 1.5707963267949;
x58.lo = -1.5707963267949; x58.up = 1.5707963267949;
x59.lo = -1.5707963267949; x59.up = 1.5707963267949;
x60.lo = -1.5707963267949; x60.up = 1.5707963267949;
x61.lo = -1.5707963267949; x61.up = 1.5707963267949;
x62.lo = -1.5707963267949; x62.up = 1.5707963267949;
x63.lo = -1.5707963267949; x63.up = 1.5707963267949;
x64.lo = -1.5707963267949; x64.up = 1.5707963267949;
x65.lo = -1.5707963267949; x65.up = 1.5707963267949;
x66.lo = -1.5707963267949; x66.up = 1.5707963267949;
x67.lo = -1.5707963267949; x67.up = 1.5707963267949;
x68.lo = -1.5707963267949; x68.up = 1.5707963267949;
x69.lo = -1.5707963267949; x69.up = 1.5707963267949;
x70.lo = -1.5707963267949; x70.up = 1.5707963267949;
x71.lo = -1.5707963267949; x71.up = 1.5707963267949;
x72.lo = -1.5707963267949; x72.up = 1.5707963267949;
x73.lo = -1.5707963267949; x73.up = 1.5707963267949;
x74.lo = -1.5707963267949; x74.up = 1.5707963267949;
x75.lo = -1.5707963267949; x75.up = 1.5707963267949;
x76.lo = -1.5707963267949; x76.up = 1.5707963267949;
x77.lo = -1.5707963267949; x77.up = 1.5707963267949;
x78.lo = -1.5707963267949; x78.up = 1.5707963267949;
x79.lo = -1.5707963267949; x79.up = 1.5707963267949;
x80.lo = -1.5707963267949; x80.up = 1.5707963267949;
x81.lo = -1.5707963267949; x81.up = 1.5707963267949;
x82.lo = -1.5707963267949; x82.up = 1.5707963267949;
x83.lo = -1.5707963267949; x83.up = 1.5707963267949;
x84.lo = -1.5707963267949; x84.up = 1.5707963267949;
x85.lo = -1.5707963267949; x85.up = 1.5707963267949;
x86.lo = -1.5707963267949; x86.up = 1.5707963267949;
x87.lo = -1.5707963267949; x87.up = 1.5707963267949;
x88.lo = -1.5707963267949; x88.up = 1.5707963267949;
x89.lo = -1.5707963267949; x89.up = 1.5707963267949;
x90.lo = -1.5707963267949; x90.up = 1.5707963267949;
x91.lo = -1.5707963267949; x91.up = 1.5707963267949;
x92.lo = -1.5707963267949; x92.up = 1.5707963267949;
x93.lo = -1.5707963267949; x93.up = 1.5707963267949;
x94.lo = -1.5707963267949; x94.up = 1.5707963267949;
x95.lo = -1.5707963267949; x95.up = 1.5707963267949;
x96.lo = -1.5707963267949; x96.up = 1.5707963267949;
x97.lo = -1.5707963267949; x97.up = 1.5707963267949;
x98.lo = -1.5707963267949; x98.up = 1.5707963267949;
x99.lo = -1.5707963267949; x99.up = 1.5707963267949;
x100.lo = -1.5707963267949; x100.up = 1.5707963267949;
x101.lo = -1.5707963267949; x101.up = 1.5707963267949;
x102.fx = 0;
x203.fx = 0;
x303.fx = 5;
x304.fx = 0;
x404.fx = 45;
x405.fx = 0;
x505.fx = 0;

* set non-default levels
x204.l = 0.1;
x205.l = 0.15;
x206.l = 0.2;
x207.l = 0.25;
x208.l = 0.3;
x209.l = 0.35;
x210.l = 0.4;
x211.l = 0.45;
x212.l = 0.5;
x213.l = 0.55;
x214.l = 0.6;
x215.l = 0.65;
x216.l = 0.7;
x217.l = 0.75;
x218.l = 0.8;
x219.l = 0.85;
x220.l = 0.9;
x221.l = 0.95;
x222.l = 1;
x223.l = 1.05;
x224.l = 1.1;
x225.l = 1.15;
x226.l = 1.2;
x227.l = 1.25;
x228.l = 1.3;
x229.l = 1.35;
x230.l = 1.4;
x231.l = 1.45;
x232.l = 1.5;
x233.l = 1.55;
x234.l = 1.6;
x235.l = 1.65;
x236.l = 1.7;
x237.l = 1.75;
x238.l = 1.8;
x239.l = 1.85;
x240.l = 1.9;
x241.l = 1.95;
x242.l = 2;
x243.l = 2.05;
x244.l = 2.1;
x245.l = 2.15;
x246.l = 2.2;
x247.l = 2.25;
x248.l = 2.3;
x249.l = 2.35;
x250.l = 2.4;
x251.l = 2.45;
x252.l = 2.5;
x253.l = 2.55;
x254.l = 2.6;
x255.l = 2.65;
x256.l = 2.7;
x257.l = 2.75;
x258.l = 2.8;
x259.l = 2.85;
x260.l = 2.9;
x261.l = 2.95;
x262.l = 3;
x263.l = 3.05;
x264.l = 3.1;
x265.l = 3.15;
x266.l = 3.2;
x267.l = 3.25;
x268.l = 3.3;
x269.l = 3.35;
x270.l = 3.4;
x271.l = 3.45;
x272.l = 3.5;
x273.l = 3.55;
x274.l = 3.6;
x275.l = 3.65;
x276.l = 3.7;
x277.l = 3.75;
x278.l = 3.8;
x279.l = 3.85;
x280.l = 3.9;
x281.l = 3.95;
x282.l = 4;
x283.l = 4.05;
x284.l = 4.1;
x285.l = 4.15;
x286.l = 4.2;
x287.l = 4.25;
x288.l = 4.3;
x289.l = 4.35;
x290.l = 4.4;
x291.l = 4.45;
x292.l = 4.5;
x293.l = 4.55;
x294.l = 4.6;
x295.l = 4.65;
x296.l = 4.7;
x297.l = 4.75;
x298.l = 4.8;
x299.l = 4.85;
x300.l = 4.9;
x301.l = 4.95;
x302.l = 5;
x305.l = 0.9;
x306.l = 1.35;
x307.l = 1.8;
x308.l = 2.25;
x309.l = 2.7;
x310.l = 3.15;
x311.l = 3.6;
x312.l = 4.05;
x313.l = 4.5;
x314.l = 4.95;
x315.l = 5.4;
x316.l = 5.85;
x317.l = 6.3;
x318.l = 6.75;
x319.l = 7.2;
x320.l = 7.65;
x321.l = 8.1;
x322.l = 8.55;
x323.l = 9;
x324.l = 9.45;
x325.l = 9.9;
x326.l = 10.35;
x327.l = 10.8;
x328.l = 11.25;
x329.l = 11.7;
x330.l = 12.15;
x331.l = 12.6;
x332.l = 13.05;
x333.l = 13.5;
x334.l = 13.95;
x335.l = 14.4;
x336.l = 14.85;
x337.l = 15.3;
x338.l = 15.75;
x339.l = 16.2;
x340.l = 16.65;
x341.l = 17.1;
x342.l = 17.55;
x343.l = 18;
x344.l = 18.45;
x345.l = 18.9;
x346.l = 19.35;
x347.l = 19.8;
x348.l = 20.25;
x349.l = 20.7;
x350.l = 21.15;
x351.l = 21.6;
x352.l = 22.05;
x353.l = 22.5;
x354.l = 22.95;
x355.l = 23.4;
x356.l = 23.85;
x357.l = 24.3;
x358.l = 24.75;
x359.l = 25.2;
x360.l = 25.65;
x361.l = 26.1;
x362.l = 26.55;
x363.l = 27;
x364.l = 27.45;
x365.l = 27.9;
x366.l = 28.35;
x367.l = 28.8;
x368.l = 29.25;
x369.l = 29.7;
x370.l = 30.15;
x371.l = 30.6;
x372.l = 31.05;
x373.l = 31.5;
x374.l = 31.95;
x375.l = 32.4;
x376.l = 32.85;
x377.l = 33.3;
x378.l = 33.75;
x379.l = 34.2;
x380.l = 34.65;
x381.l = 35.1;
x382.l = 35.55;
x383.l = 36;
x384.l = 36.45;
x385.l = 36.9;
x386.l = 37.35;
x387.l = 37.8;
x388.l = 38.25;
x389.l = 38.7;
x390.l = 39.15;
x391.l = 39.6;
x392.l = 40.05;
x393.l = 40.5;
x394.l = 40.95;
x395.l = 41.4;
x396.l = 41.85;
x397.l = 42.3;
x398.l = 42.75;
x399.l = 43.2;
x400.l = 43.65;
x401.l = 44.1;
x402.l = 44.55;
x403.l = 45;
x507.l = 0.01;

Model m / all /;

m.limrow=0; m.limcol=0;
m.tolproj=0.0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

$if not set NLP $set NLP NLP
Solve m using %NLP% minimizing objvar;
