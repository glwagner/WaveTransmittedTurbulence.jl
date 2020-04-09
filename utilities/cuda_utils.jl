using CUDAnative, CUDAdrv

function select_device!(ndev)
    @show dev = CuDevice(ndev)
    CUDAnative.device!(dev)
    return nothing
end
