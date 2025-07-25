using Plots

file = ARGS[1]

# parse csv and sum lanthanides for each x for each timestep, producing a histogram for every n
function sumLanthanides(gridsize)
    lncsv = open("$file", "r")
    t = 1
    sums = zeros(gridsize)

    # Import ln.csv as array
    for line in readlines(lncsv)
        if (length(line) > 0)
            if (length(line) < gridsize / 2)
                # Time
                t = parse(Float64, line)
                if (mod(t,300) == 0)
                    display(t)
                    println(sums)
                    bar(1:first(size(sums)), sums, legend=false, ylimits=(0,6000))
                    savefig("$t.png")
                end
                sums = zeros(gridsize)
            else
                # Row
                j = 1
                for v in split(line, ",")
                    sums[j] += parse(Float64, v)
                    j+=1
                end
            end
        end
    end
end

sumLanthanides(100)
