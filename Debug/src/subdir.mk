################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/pic2dCluster.cu 

CU_DEPS += \
./src/pic2dCluster.d 

OBJS += \
./src/pic2dCluster.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-6.5/bin/nvcc -I/usr/local/lib -G -g -O0 -gencode arch=compute_35,code=sm_35 --target-cpu-architecture x86 -m64 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-6.5/bin/nvcc -I/usr/local/lib -G -g -O0 --compile --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 --target-cpu-architecture x86 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


