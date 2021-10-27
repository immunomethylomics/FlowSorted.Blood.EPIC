.onLoad <- function(libname = find.package("FlowSorted.Blood.EPIC"), 
                    pkgname = "FlowSorted.Blood.EPIC"){
    
    
    if(getRversion() >= "3.4.0") 
        utils::globalVariables(c("RGsetTargets"))
    invisible()
    
    titles <- strwrap("FlowSorted.Blood.EPIC: Illumina Human Methylation data  
                        from EPIC on immunomagnetic sorted adult blood cell 
                        populations",
                        width = 128)
    createHubAccessors(pkgname, titles)
    
}