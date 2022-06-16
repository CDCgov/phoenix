/*
========================================================================================
    Processes
========================================================================================
*/



process database_check {
    db_ch = Channel.fromPath(${params.databases}, checkIfExists: true )

    input:
    path(db_path)

    output:

    script:
    """
        database_checker.sh ${db_path}
    """
}

workflow database_check {
  database_check(params.databases)
}
