//
// This workflow will handles taking in either a samplesheet or directory and creates correct channels.
//

workflow ADD_EMPTY_CHANNEL {
    take:
        incoming_ch  // incoming channel

    main:
        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        output_ch = incoming_ch.map{ add_empty_ch(it) } // create meta for sample

    emit:
        output_ch = output_ch       // channel: [ meta, [ scaffolds_file ] ]

}

def add_empty_ch(input_ch) {
    meta_id = input_ch[0]
    output_array = [ meta_id, input_ch[1], input_ch[2], []]
    return output_array
}