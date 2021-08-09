// Load pipeline configurations such as reference genome, locations of tools
strling='strling'
outliers_py='strling-outliers.py'
load "pipeline_config.groovy"

if(args.any { it.endsWith('.cram') })
    input_type = 'cram'
else
    input_type='bam'

def get_fname(path) {
    def x = path.split("/")[-1]
    return(x)
}

str_index = {
    def str_ref = get_fname(REF) + ".str"
    produce(str_ref) {
        exec """
            $strling index
                $REF
                -g $str_ref
        ""","strling"
    }
}

str_extract = {
    def sample = branch.name
    def str_ref = get_fname(REF) + ".str"
    produce(sample + ".str.bin") {
        exec """
            $strling extract
                -f $REF
                -g $str_ref
                ${input[input_type]}
                $output.bin
        ""","strling"
    }
}

str_merge = {
    from("*.bin") produce("strling-bounds.txt") {
        exec """
            $strling merge
                -f $REF
                $inputs.bin
        ""","merge"
    }
}

str_call_individual = {
    def sample = branch.name.prefix
    from (sample + "*str.bin", sample + "*" + input_type) produce(
            sample + "-bounds.txt", sample + "-unplaced.txt",
            sample + "-genotype.txt") {
        exec """
            $strling call
                -f $REF
                -o $sample
                ${input[input_type]}
                $input.bin
        ""","strling"
    }
}

str_call_joint = {
    def sample = branch.name.prefix
    from (sample + "*str.bin", sample + "*" + input_type, "strling-bounds.txt") produce(
            sample + "-bounds.txt", sample + "-unplaced.txt",
            sample + "-genotype.txt") {
        exec """
            $strling call
                -f $REF
                -b $input.txt
                -o $sample
                ${input[input_type]}
                $input.bin
        ""","joint_call"
    }
}

outliers = {
    from("*-genotype.txt", "*-unplaced.txt") produce("STRs.tsv") {
        exec """
            $outliers_py
                --genotypes *-genotype.txt
                --unplaced *-unplaced.txt
                --emit control-file.tsv
        """
    }
}
