MODULE array_bound_module

  use parallel_module, only: id_grid,ir_grid,id_band,ir_band &
                            ,id_bzsm,ir_bzsm,id_spin,ir_spin &
                            ,myrank_g,myrank_b,myrank_k,myrank_s

  implicit none

  PRIVATE
  PUBLIC :: set_array_bound

  integer,PUBLIC :: ML ,ML_0 ,ML_1
  integer,PUBLIC :: MB ,MB_0 ,MB_1
  integer,PUBLIC :: MBZ,MBZ_0,MBZ_1
  integer,PUBLIC :: MSP,MSP_0,MSP_1

CONTAINS

  SUBROUTINE set_array_bound
    ML_0  = id_grid(myrank_g)+1
    ML_1  = id_grid(myrank_g)+ir_grid(myrank_g)
    MB_0  = id_band(myrank_b)+1
    MB_1  = id_band(myrank_b)+ir_band(myrank_b)
    MBZ_0 = id_bzsm(myrank_k)+1
    MBZ_1 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)
    MSP_0 = id_spin(myrank_s)+1
    MSP_1 = id_spin(myrank_s)+ir_spin(myrank_s)
    ML  = sum(ir_grid)
    MB  = sum(ir_band)
    MBZ = sum(ir_bzsm)
    MSP = sum(ir_spin)
  END SUBROUTINE set_array_bound

END MODULE array_bound_module
