#object = ret$infect
#molecules = prot
#selected.molecules=character(0)
#additional.variables=NULL
#parameters=ewaff.report.parameters()
#verbose=T
msg <- function(..., verbose=T) {
    if (verbose) {
        x <- paste(list(...))
        name <- sys.call(sys.parent(1))[[1]]
        cat(paste("[", name, "]", sep=""), date(), x, "\n")
    }
}

#' Generate CpG associations report
#'
#' Generate HTML file that summarises CpG site associations. 
#'
#' @param  object Output from \code{\link{ewaff.summary}}.
#' @param  output.file Default = "ewas-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "Illumina methylation data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' 
#' @export
ewaff.report <- function(object,
                         output.file = "report.html",
                         author = "Analyst",
                         study = "Illumina methylation data",
                         ...) {
    stopifnot(is.summary.object(object))
    
    ewaff:::msg("Writing report as html file to", output.file)
    report.path <- system.file("report", package="ewaff")
    require(knitr)
    require(Cairo)
    require(gridExtra)

    options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
    
    opts <- opts_chunk$get()
    on.exit(opts_chunk$set(opts))
    opts_chunk$set(warning=FALSE, echo=FALSE, message=FALSE, results="asis",
                   fig.width=6, fig.height=6, dev="CairoPNG")
    knit.report(file.path(report.path, "report.rmd"),output.file, ...)
}


#' Summarize results.
#'
#' Generates variable and covariate summary tables, QQ plots,
#' Manhattan plots, a list of associations, plots of the strongest
#' associations and plots of selected CpG sites.
#'
#' @param object Object returned by \code{\link{ewaff.sites}}.
#' @param molecules Matrix of molecules levels used in the analysis.
#' @param selected.molecules Vector of CpG site names to plot (Default: character(0)).
#' @param additional.variables A data frame of variables to show associations with the variable of interest.
#' @param parameters Default = ewaff.report.parameters(). List of parameter values.
#' See \code{\link{ewaff.report.parameters}()}.
#' @export
#' @return List
#' 
#' @export
prot.summary <- function(object, molecules,
                          selected.molecules=character(0),
                          additional.variables=NULL,
                          parameters=ewaff.report.parameters(),
                          verbose=T) {
    
    stopifnot(ewaff:::is.sites.object(object))

    stopifnot(is.null(additional.variables) || nrow(additional.variables) == ncol(molecules))
    
    if (ncol(molecules) > length(object$sample.idx)) {
        molecules <- molecules[,object$sample.idx]
        if (!is.null(additional.variables))
            additional.variables <- additional.variables[object$sample.idx,,drop=F]
    }
    
    p.values <- object$table$p.value
    p.adjusted <- object$table$p.adjust
    estimates <- if (!is.null(object$table$estimate)){ 
                    object$table$estimate 
                } else object$table$f
    
    stopifnot(length(p.values) == nrow(molecules))
    stopifnot(parameters$max.plots < length(p.values))    
    stopifnot(all(selected.molecules %in% rownames(molecules)))
    
    parameters$practical.threshold <- p.values[order(p.values)[parameters$max.plots+1]]

    if (is.na(parameters$sig.threshold)) 
        parameters$sig.threshold <- 0.05/nrow(molecules)

    sig.idx <- which(p.values < parameters$sig.threshold)
    practical.idx <- which(p.values < parameters$practical.threshold)
    selected.idx <- match(selected.molecules, rownames(molecules))    

    mol.idx<- union(practical.idx, selected.idx)

    mol.stats <- data.frame(estimate=estimates[mol.idx],
                            p.value=p.values[mol.idx],
                            p.adjust=p.adjusted[mol.idx])

    rownames(mol.stats) <- rownames(molecules)[mol.idx]
    
    ewaff:::msg("QQ plots", verbose=verbose)
    qq.plot <- ewaff.qq.plot(p.values=p.values,
                              sig.threshold=parameters$sig.threshold,
                              lambda.method=parameters$qq.inflation.method)

    msg("volcano plots", verbose=verbose)
    volcano.plot <- data.frame(p.values=p.values,
                                estimates=estimates) |>
                        ggplot(aes(x = estimates, y = -log10(p.values))) +
                        geom_point() +
                        geom_hline(yintercept = -log10(parameters$sig.threshold), 
                            linetype='dashed')

    plot.sites <- rownames(molecules)[union(practical.idx, selected.idx)]
    msg("CpG site plots:", length(plot.sites), verbose=verbose)
    variable.of.interest <- ifelse(object$independent.variable == "methylation",
                                   object$dependent.variable, object$independent.variable)

    if (!is.null(object$table$estimate)){ 
        mol.plots <- sapply(plot.sites, function(cpg) {
            msg("Plotting", cpg, verbose=verbose)
            ewaff:::ewaff.cpg.plot(variable.of.interest, as.data.frame(object$design), molecules[cpg,], cpg)
        }, simplify=F)
    } else mol.plots <- NULL

    sample.characteristics <- NULL
    covariate.associations <- NULL
    additional.associations <- NULL
    #if (object$method != "coxph") {
    #    msg("Sample characteristics", verbose=verbose)
    #    colnames.to.keep <- colnames(object$design)[-1]
    #    data <- as.data.frame(object$design[,-1], drop = F)
    #    colnames(data) <- colnames.to.keep
    #    sample.characteristics <- ewaff.sample.characteristics(variable.of.interest, data)
    #    covariate.associations <- ewaff.covariate.associations(variable.of.interest, data)
    #    if (!is.null(additional.variables)) {
    #        additional.variables <- cbind(data[,variable.of.interest,drop=F], additional.variables)
    #        additional.associations <- ewaff.covariate.associations(variable.of.interest, additional.variables)
    #    }
    #}
    ##parameters$winsorize.pct <- object$winsorize.pct
    
    list(class="ewaff.summary",
         parameters=parameters,
         qq.plot=qq.plot,
         volcano.plot=volcano.plot,
         mol.stats=mol.stats,
         mol.plots=mol.plots,
         practical.sites=rownames(molecules)[practical.idx],
         significant.sites=rownames(molecules)[sig.idx],
         selected.sites=rownames(molecules)[selected.idx],         
         sample.characteristics=sample.characteristics,
         covariate.associations=covariate.associations,
         additional.associations=additional.associations)
}

is.summary.object <- function(object)
    is.list(object) && "class" %in% names(object) && object$class == "ewaff.summary"
