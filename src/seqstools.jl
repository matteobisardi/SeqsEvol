export num2letter, vec2string

"""
	num2letter(i)

	Converts the amino acid "i" into its corresponding letter
	in the mapping "ACD...Y-"  <---> "1 2 3 ... 20 21".
"""

let alphabet = ["A", "C", "D", "E", "F", "G", "H", "I",  "K", "L",  "M",  "N", "P",  "Q",  "R",
"S",  "T", "V",  "W",  "Y"]
    global num2letter
    function num2letter(i :: Integer)
        1 <= i <= 20 && return alphabet[i]
        return "-"
    end
end


"""
	vec2string(v)

	Converts the the vector "v" of amino acids in number format
	to its equivalent in letters following the mapping "ACD...Y-"  <---> "1 2 3 ... 20 21".
"""

function vec2string(v)
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end
