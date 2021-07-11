using Logging
using LoggingExtras
path = "C:/Users/riw/tubCloud/Uni/Market_Tool/pomato_studies/data_temp/julia_files/data/cwe_2030/log.txt"

function logshit()
    for i in 1:10
        @info(i, "logogo")
    end
end
io = open(path, "w+")
logger = FileLogger(io)
global_logger(logger)
@info("a context specific log message")
logshit()

a = Dict("a" => 1, "b" => 2)
for k in keys(a)
    println(k, a[k])
end