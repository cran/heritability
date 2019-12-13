repeatability <-
function (data.vector, geno.vector, line.repeatability = FALSE,
    covariates.frame = data.frame())
{

# data.vector=d$y; geno.vector=d$genotype; covariates.frame = data.frame(d$x); line.repeatability = FALSE

    stopifnot(length(geno.vector) == length(data.vector))
    her.frame <- data.frame(dat = data.vector, geno = geno.vector)
    her.frame <- her.frame[!is.na(her.frame$dat), ]
    number.of.covariates <- 0
    if (nrow(covariates.frame) > 0) {
        stopifnot(nrow(covariates.frame) == length(data.vector))
        cov.names <- names(covariates.frame)
        number.of.covariates <- length(cov.names)
        covariates.frame <- as.data.frame(covariates.frame[!is.na(data.vector),
            ])
        names(covariates.frame) <- cov.names
        her.frame <- cbind(her.frame, covariates.frame)
        her.frame <- her.frame[apply(covariates.frame, 1, function(x) {
            sum(is.na(x))
        }) == 0, ]
    }
    her.frame$geno <- factor(her.frame$geno)
    n.rep.vector <- as.integer(table(her.frame$geno))
    n.geno <- length(n.rep.vector)
    average.number.of.replicates <- (sum(n.rep.vector) - sum(n.rep.vector^2)/sum(n.rep.vector))/(length(n.rep.vector) -
        1)
    if (max(n.rep.vector) == 1) {
        return(list(repeatability = NA, gen.variance = NA, res.variance = NA))
    }
    else {
        if (nrow(covariates.frame) > 0) {
            av <- anova(lm(as.formula(paste("dat~", paste(cov.names,
                collapse = "+"), "+geno")), data = her.frame))
        }
        else {
            av <- anova(lm(dat ~ geno, data = her.frame))
        }
        gen.variance <- (av[[3]][1 + number.of.covariates] - av[[3]][2 + number.of.covariates])/average.number.of.replicates
        if (gen.variance < 0) {
            gen.variance <- 0
        }
        if (line.repeatability) {
            res.variance <- av[[3]][2 + number.of.covariates]/average.number.of.replicates
        }
        else {
            res.variance <- av[[3]][2 + number.of.covariates]
        }
        if (!line.repeatability) {
            F.ratio <- av[[3]][1 + number.of.covariates]/av[[3]][2 + number.of.covariates]
            df.1 <- n.geno - 1
            df.2 <- av[[1]][2 + number.of.covariates]
            F.L <- qf(0.025, df1 = df.1, df2 = df.2, lower.tail = TRUE)
            F.U <- qf(0.975, df1 = df.1, df2 = df.2, lower.tail = TRUE)
            conf.int.left <- (F.ratio/F.U - 1)/(F.ratio/F.U +
                average.number.of.replicates - 1)
            conf.int.right <- (F.ratio/F.L - 1)/(F.ratio/F.L +
                average.number.of.replicates - 1)
            if (conf.int.left < 0) {
                conf.int.left <- 0
            }
            if (conf.int.right > 1) {
                conf.int.right <- 1
            }
            if (conf.int.left > 1) {
                conf.int.left <- 1
            }
            if (conf.int.right < 0) {
                conf.int.right <- 0
            }
        }
        else {
            conf.int.left <- NA
            conf.int.right <- NA
        }
        return(list(repeatability = gen.variance/(gen.variance +
            res.variance), gen.variance = gen.variance, res.variance = res.variance,
            line.repeatability = line.repeatability, average.number.of.replicates = average.number.of.replicates,
            conf.int = c(conf.int.left, conf.int.right)))
    }
}