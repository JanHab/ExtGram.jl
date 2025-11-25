# simple_script.jl

# Check if any arguments were passed
if isempty(ARGS)
    println("Please provide some arguments!")
    exit(1)
end

# Access arguments by index
first_arg = ARGS[1]

# Remember: All arguments come in as Strings.
# You must parse them if you need numbers.
if length(ARGS) >= 2
    second_arg = parse(Int, ARGS[2])
    result = second_arg * 2
else
    second_arg = "N/A"
    result = 0
end

println("Argument 1 (String): $first_arg")
println("Argument 2 (Number): $second_arg")
println("Calculation Result: $result")