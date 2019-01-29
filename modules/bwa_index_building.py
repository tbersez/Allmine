#bwa index building

rule bwa_index:
    input:
        genome = config["REF"] + config["GENOME"]
    output:
        protected(config["REF"] + config["GENOME"] + ".amb"),
        protected(config["REF"] + config["GENOME"] + ".ann"),
        protected(config["REF"] + config["GENOME"] + ".bwt"),
        protected(config["REF"] + config["GENOME"] + ".pac"),
        protected(config["REF"] + config["GENOME"] + ".sa")
    message: "Building BWA index for reference genome {input.genome}\n"
    threads: config["THREADS"]
    shell:
        """
        /home/aa/anaconda3/bin/bwa index\
        {input.genome}
        """
