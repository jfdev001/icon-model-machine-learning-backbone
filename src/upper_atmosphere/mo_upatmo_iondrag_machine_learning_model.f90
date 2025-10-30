! @file mo_upatmo_iondrag_machine_learning_model.f90
! 
! @brief Define a "hello world" FTorch routine.
!
! @references 
! * https://github.com/Cambridge-ICCS/FTorch/blob/ef44cf6d70edef38003dec41ee7a1b496922a1a5/examples/1_Tensor/tensor_manipulation.f90
MODULE mo_upatmo_iondrag_machine_learning_model
  USE ftorch, ONLY: assignment(=), operator(+), torch_kCPU, torch_kFloat32, &
    torch_tensor, torch_tensor_delete, torch_tensor_ones, &
    torch_tensor_from_array

  USE, INTRINSIC :: iso_c_binding, only: c_int64_t

  ! Import the real32 TYPE for 32-bit floating point numbers
  USE, INTRINSIC :: iso_fortran_env, only: sp => real32

  IMPLICIT NONE

  PUBLIC :: hello_ftorch

CONTAINS 

  SUBROUTINE hello_ftorch(dst)
      INTEGER, INTENT(IN) :: dst 

      ! Set working precision for reals to be 32-bit floats
      INTEGER, PARAMETER :: wp = sp

      ! Define some tensors
      TYPE(torch_tensor) :: a, b, c

      ! Variables for constructing tensors with torch_tensor_ones
      INTEGER, PARAMETER :: ndims = 2
      INTEGER(c_int64_t), DIMENSION(2), PARAMETER :: tensor_shape = [2, 3]

      ! Variables for constructing tensors with torch_tensor_from_array
      REAL(wp), DIMENSION(2,3), TARGET :: in_data, out_data

      ! Create tensors, perform simple operations, then output result
      WRITE (dst, '(3a)') "mo_iondrag_machine_learning_model: hello ftorch"

      CALL torch_tensor_ones(a, ndims, tensor_shape, torch_kFloat32, torch_kCPU)
      in_data(:,:) = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [2,3])
      CALL torch_tensor_from_array(b, in_data, torch_kCPU)
      CALL torch_tensor_from_array(c, out_data, torch_kCPU)

      c = a + b
      WRITE(dst, '(3a)') "mo_iondrag_machine_learning_model: sum of tensors"
      WRITE(dst, *) out_data
      WRITE(dst, '(3a)') "mo_iondrag_machine_learning_model: expected output: reshape([2, 3, 4, 5, 6, 7], [2,3])"

      CALL torch_tensor_delete(a)
      CALL torch_tensor_delete(b)
      CALL torch_tensor_delete(c)
  END SUBROUTINE hello_ftorch

END MODULE mo_upatmo_iondrag_machine_learning_model
