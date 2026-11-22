"""String rewriting systems from Yolcu-Aaronson-Heule (arXiv 2105.14697v3).

ASCII symbol names:
    '<'  = the left marker (paper: ◁), functional view <(x) = 1
    '>'  = the right marker (paper: ⊿), functional view >(x) = x
    '*'  = the odd right marker (paper: •), *(x) = 2x + 1
    'f'  = binary 0,   f(x) = 2x
    't'  = binary 1,   t(x) = 2x + 1
    '0','1','2' = ternary digits, d(x) = 3x + d
A string s1 s2 ... sn denotes sn(...(s2(s1(.)))) with s1 = '<' the constant 1.

For the 5x+1 control system the digits '0'..'4' are base-5 digits d(x) = 5x + d.
"""

# ---------------------------------------------------------------- auxiliary X = A u B
A = [("f0", "0f"), ("f1", "0t"), ("f2", "1f"),
     ("t0", "1t"), ("t1", "2f"), ("t2", "2t")]
B = [("<0", "<t"), ("<1", "<ff"), ("<2", "<ft")]
X = A + B

# ---------------------------------------------------------------- dynamic rule sets
D_T = [("f>", ">"), ("t>", "2>")]                      # T(n) = n/2 | (3n+1)/2
D_W = [("f>", "0>")]                                   # W': 2x -> 3x  (Thm 4.1)
D_F = [("1>", ">"), ("0f>", "0>"), ("1f>", "1>"),      # Farkas variant F (Sec 4.2)
       ("1t>", "12>"), ("2t>", "22>")]
D_S = [("ff*", "0*"), ("tf*", "*"), ("t*", "2*")]      # accelerated S (Sec 4.4)
D_EVEN4 = [("ff>", "f>"), ("tf>", "t>"), ("t>", "2f>")]          # C, mod 4
D_ODD4 = [("f>", ">"), ("ft>", "1f>"), ("tt>", "2t>")]           # T, mod 4
D_EVEN8 = [("fff>", "ff>"), ("ftf>", "ft>"), ("tff>", "tf>"),    # C, mod 8
           ("ttf>", "tt>"), ("t>", "2f>")]
D_ODD8 = [("f>", ">"), ("fft>", "0tf>"), ("ftt>", "1ft>"),       # T, mod 8
          ("tft>", "2ff>"), ("ttt>", "2tt>")]

# ---------------------------------------------------------------- 5x+1 analogue (binary/quinary)
# swap rules: (binary b)(quinary q) = 10x + 5b + q  ->  (quinary q')(binary b') = 10x + 2q' + b'
A5 = []
for b, bs in ((0, "f"), (1, "t")):
    for q in range(5):
        v = 5 * b + q
        A5.append((bs + str(q), str(v // 2) + ("t" if v % 2 else "f")))
# <q has value 5 + q; write 5+q in binary without the leading 1
B5 = [("<" + str(q), "<" + bin(5 + q)[3:].replace("0", "f").replace("1", "t")) for q in range(5)]
X5 = A5 + B5
D_T5 = [("f>", ">"), ("t>", "3>")]      # (5(2x+1)+1)/2 = 5x+3


def remove(rules, *lhs):
    out = [r for r in rules if r[0] not in lhs]
    assert len(out) == len(rules) - len(lhs), lhs
    return out


SYSTEMS = {
    "T": D_T + X,
    "Wp": D_W + X,
    "F": D_F + X,
    "S": D_S + X,
    "S1": remove(D_S + X, "ff*"),       # open in YAH (Prop 4.5.1)
    "S2": remove(D_S + X, "tf*"),       # open in YAH (Prop 4.5.2)
    "S3": remove(D_S + X, "t*"),        # proved in YAH (Prop 4.6)
    "EVEN8": D_EVEN8 + X,
    "ODD8": D_ODD8 + X,
    "T5": D_T5 + X5,                    # 5x+1 analogue: NOT terminating (cycle 13 -> 33 -> ...)
}
for lhs, _ in D_T + X:
    SYSTEMS["T-" + lhs] = remove(D_T + X, lhs)          # Table 1: all 11 proved by YAH
for lhs, _ in D_EVEN4:
    SYSTEMS["EVEN4-" + lhs] = remove(D_EVEN4 + X, lhs)  # Thm 4.9 (proved)
for lhs, _ in D_ODD4:
    SYSTEMS["ODD4-" + lhs] = remove(D_ODD4 + X, lhs)    # Thm 4.9 (proved)
for lhs, _ in D_EVEN8:
    SYSTEMS["EVEN8-" + lhs] = remove(D_EVEN8 + X, lhs)  # 0 mod 8 (fff>) open; 2,4,6 proved
for lhs, _ in D_ODD8:
    SYSTEMS["ODD8-" + lhs] = remove(D_ODD8 + X, lhs)    # 1,5,7 mod 8 open; 3 proved
SYSTEMS["ODD8-tft>-ttt>"] = remove(D_ODD8 + X, "tft>", "ttt>")  # Conj 4.14 (open)
SYSTEMS["ODD8-fft>-tft>"] = remove(D_ODD8 + X, "fft>", "tft>")  # {1,5}: proved by YAH
SYSTEMS["ODD8-fft>-ttt>"] = remove(D_ODD8 + X, "fft>", "ttt>")  # {1,7}: proved by YAH
SYSTEMS["T5-t>"] = remove(D_T5 + X5, "t>")                       # trivially terminating

# ---- accelerated encodings of Conjecture 4.14 (not tried as SRSs in YAH, as far as the paper says)
# H(n) = 3n/4 (n = 0 mod 4), (9n+1)/8 (n = 7 mod 8), undefined otherwise (YAH eq. (10), Prop 4.15)
D_H = [("ff>", "0>"), ("ttt>", "22>")]
SYSTEMS["H"] = D_H + X
# G(n) = (3n+1)/4 (n = 1 mod 32), (9n+7)/8 (n = 57 mod 64)  (YAH proof of Prop 4.15)
D_G = [("fffft>", "0fft>"), ("tttfft>", "22fft>")]
SYSTEMS["G"] = D_G + X
_EXTRA_STATUS = {"H": "O(=Conj 4.14)", "G": "O(=Conj 4.14)"}

# ---- 5x+1 analogue of the mod-8 odd system (T5(n) = n/2 | (5n+1)/2), binary/quinary
# rhs for 8x+r: 20x + (5r+1)/2 = 5(4x) + (5r+1)/2 - 5k ... written as "ff" + quinary digit
# when (5r+1)/2 < 5, else "fft" etc. computed below by value.
def _t5_rhs(r):
    c = (5 * r + 1) // 2           # T5(8x+r) = 20x + c, c in {3, 8, 13, 18}
    q, d = divmod(c, 5)            # 20x + c = 5(4x + q) + d
    pre = {0: "ff", 1: "ft", 2: "tf", 3: "tt"}[q]   # 4x + q
    return pre + str(d) + ">"
D_ODD8_5 = [("f>", ">")] + [(bits + ">", _t5_rhs(r)) for bits, r in
                            (("fft", 1), ("ftt", 3), ("tft", 5), ("ttt", 7))]
SYSTEMS["ODD8_5"] = D_ODD8_5 + X5
SYSTEMS["ODD8_5-tft>-ttt>"] = remove(D_ODD8_5 + X5, "tft>", "ttt>")  # FALSE analogue: cycle 17
_EXTRA_STATUS.update({"ODD8_5": "N", "ODD8_5-tft>-ttt>": "N"})

DIGIT_BASE = {"T5": 5, "T5-t>": 5, "ODD8_5": 5, "ODD8_5-tft>-ttt>": 5}

# Status as reported in the paper (P = proved terminating there, O = open, N = nonterminating)
YAH_STATUS = {"T": "O", "Wp": "P", "F": "P", "S": "O", "S1": "O", "S2": "O", "S3": "P",
              "T5": "N", "ODD8-tft>-ttt>": "O", "ODD8-fft>-tft>": "P", "ODD8-fft>-ttt>": "P",
              "EVEN8-fff>": "O", "EVEN8-ftf>": "P", "EVEN8-tff>": "P", "EVEN8-ttf>": "P",
              "ODD8-fft>": "O", "ODD8-ftt>": "P", "ODD8-tft>": "O", "ODD8-ttt>": "O"}
YAH_STATUS.update(_EXTRA_STATUS)
for k in SYSTEMS:
    if k.startswith(("T-", "EVEN4-", "ODD4-")):
        YAH_STATUS[k] = "P"



def head_markers(rules):
    """Rules whose lhs starts with '<' (forward top rules) and rules whose lhs ends
    with an end marker (reverse top rules)."""
    fwd = [r for r in rules if r[0].startswith("<")]
    rev = [r for r in rules if r[0][-1] in ">*"]
    return fwd, rev


def reverse(rules):
    return [(l[::-1], r[::-1]) for l, r in rules]


# ---------------------------------------------------------------- functional semantics
def apply_sym(c, v, base3=3):
    if c == "f":
        return 2 * v
    if c == "t":
        return 2 * v + 1
    if c == ">":
        return v
    if c == "*":
        return 2 * v + 1
    if c.isdigit():
        return base3 * v + int(c)
    raise ValueError(c)


def value(s, base3=3):
    assert s[0] == "<"
    v = 1
    for c in s[1:]:
        v = apply_sym(c, v, base3)
    return v


def gamma(s, x, base3=3):
    """Value of a marker-free (or '<'-free) word applied to x."""
    v = x
    for c in s:
        v = 1 if c == "<" else apply_sym(c, v, base3)
    return v
