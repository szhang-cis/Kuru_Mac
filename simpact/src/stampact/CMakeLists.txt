# stampact library
set (stampact_src ${src_stampact}/bonds_db.f90
                  ${src_stampact}/cont_db.f90
                  ${src_stampact}/data_db.f90
                  ${src_stampact}/dflib.f90
                  ${src_stampact}/flc_db.f90
                  ${src_stampact}/param_db.f90)

add_library (${stampactlib} STATIC ${stampact_src})
