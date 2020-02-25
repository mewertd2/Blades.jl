# Copyright 2020 Digital Domain 3.0
#
# Licensed under the Apache License, Version 2.0 (the "Apache License")
# with the following modification; you may not use this file except in
# compliance with the Apache License and the following modification to it:
# Section 6. Trademarks. is deleted and replaced with:
#
# 6. Trademarks. This License does not grant permission to use the trade
#    names, trademarks, service marks, or product names of the Licensor
#    and its affiliates, except as required to comply with Section 4(c) of
#    the License and to reproduce the content of the NOTICE file.
#
# You may obtain a copy of the Apache License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Apache License with the above modification is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the Apache License for the specific
# language governing permissions and limitations under the Apache License.

using Test
using Blades
using LinearAlgebra
using Combinatorics

module g3
  using Blades
  @generate_basis("+++")
end

module pg2
  using Blades
  @generate_basis("0++")
end

module g4
  using Blades
  @generate_basis("++++", false, true)
end

module pg3
  using Blades
  @generate_basis("0+++")
end

module g2
  using Blades
  @generate_basis("++")
end

module G5
  using Blades
  @generate_basis("+++++",false,true,true)
end

using .pg3
using .g4
using .g3
using .pg2
using .g2
using .G5

partial(f,x...) = (y...)->f(x...,y...)


"""
We can complete the definition of the geometric product.  
All Blades are simple k-vectors.  So only one of a∧b or a⋅b is ever nonzero. 
By defining addition of a Blade with 0 we can complete the Geometric Algebra of Blades.
u*v == u∧v + u⋅v
Not defined in the Blades package since it is defined in Multivectors
"""
Base.:(+)(a::B, s::Z) where {B<:Blade,Z<:Val{0.0}} = a
Base.:(+)(a::B, s::T) where {B<:Blade,T<:Real} = a + Val(Float64(s)) 
Base.:(+)(s::T, a::B) where {B<:Blade,T<:Real} = a + s

@testset "PGA3D" begin
  e₁, e₂, e₃, e₄, e₁₂, e₁₃, e₁₄, e₂₃, e₂₄, e₃₄, e₁₂₃, e₁₂₄, e₁₃₄, e₂₃₄, e₁₂₃₄ = alle(pg3,4)
   
  a = e₁(1.0); b = e₂(2.0); c = e₃(3.0); d = e₄(4.0)

  @test e₁(1.0)*e₁(1.0) == 0.0
  @test e₂(1)*e₂(1) == 1
  @test a*b*a == 0.0
  @test a*b == -b*a
  @test a*b*c == b*c*a 
  @test b*c*b == -b.x*b.x*c
  @test a*b == e₁₂(a.x*b.x)
  @test reverse(a*b*c*d) == d*c*b*a == reverse((a*b)*(c*d)) == (d*c*b)*a
  @test e₃(2.0)*e₂₃₄(3.0) == -e₂₄(6.0) == e₂₄(-6.0) 
  @test grade(1) == 0
  @test grade(a*b*a) == 0
  @test grade(a*b*c) == 3

  B = b*c

  @test a*B == a∧B
  @test b*B == b⋅B
  @test d*B == d∧B
  @test a*b == a∧b + a⋅b
  @test c*B == c∧B + c⋅B
  @test B*a == B∧a + B⋅a
  @test c*c == c∧c + c⋅c

  @test 1.0∧a == a
  @test 1.0⋅a == a

  g₁, g₂, g₃, g₄, g₁₂, g₁₃, g₁₄, g₂₃, g₂₄, g₃₄, g₁₂₃, g₁₂₄, g₁₃₄, g₂₃₄, g₁₂₃₄ = alle(g4,4)

  @test dual(dual(1g₁)) == ⟂(⟂(1g₁))
  @test subspace(dual(dual(1e₁))) == subspace(⟂(⟂(1g₁)))
  @test subspace(Blades.nonmetric_dual(Blades.nonmetric_dual(3e₁))) == subspace(⟂(⟂(3g₁)))
  @test sign(dual(dual(1pg2.e₁))) == sign(⟂(⟂(1g3.e₁)))
  @test sign(Blades.nonmetric_dual(Blades.nonmetric_dual(3pg2.e₁))) == sign(⟂(⟂(3g3.e₁)))
  @test sign(dual(3pg2.e₁))== sign(⟂(3g3.e₁))
  @test sign(Blades.nonmetric_dual(3pg2.e₁)) == sign(⟂(3g3.e₁))

  # convenience operators on type
  @test e₁∧e₂ == e₁₂ == e₁*e₂
  @test e₂*e₁ == e₂∧e₁ == Nothing
  @test 2(e₁∧e₂) == e₁₂(2)
end

@testset "G2" begin

  e₁, e₂, e₁₂ = alle(g2,2)
  eᵢ  = one.(alle( g2, 2)[1:2])
  #More generally, on any flat space we have
  #⋆(dxi1 ∧dxi2 ∧···∧dxik) = dxik+1 ∧dxik+2 ∧···∧dxin, where (i1,i2,...,in) is any even permutation of (1,2,...,n)
  
  for p in Iterators.filter(p->parity(p)==0,permutations(1:2))
    for s in 1:1
      @test ⋆(reduce(∧,eᵢ[p[1:s]])) == reduce(∧,eᵢ[p[s+1:end]])
    end
  end
  @test ⋆(1e₁) == 1e₂
  @test ⋆(1e₂) == -1e₁
  ot = -1.0e₁₂
  @test ot∧⋆(ot) == (ot⋅ot)pseudoscalar(ot)
  @test ⋆(2.0e₁₂) == -2.0
  ot = 2.0e₁₂
  @test ot∧⋆(ot) == (ot⋅ot)pseudoscalar(ot)
end

#==
@testset "DEC course notes G5" begin
  eᵢ  = one.(alle( G5, 5)[1:5])
  #More generally, on any flat space we have
  #⋆(dxi1 ∧dxi2 ∧···∧dxik) = dxik+1 ∧dxik+2 ∧···∧dxin, where (i1,i2,...,in) is any even permutation of (1,2,...,n)
  
  for p in Iterators.filter(p->parity(p)==0,permutations(1:5))
    for s in 1:4
      @test ⋆(reduce(∧,eᵢ[p[1:s]])) == reduce(∧,eᵢ[p[s+1:end]])
      w = reduce(∧,eᵢ[p[1:s]])
    end
  end
end
==#

@testset "hodge star" begin
  rb = alle(G5, 5)[1:5]
  𝑖 = mapreduce(eᵢ->1.0*eᵢ, ∧, rb[1:4])
  for i in 1:100
    A = rand()rand(rb[1:3])∧rand()rand(rb[1:3])
    B = rand()rand(rb[1:3])∧rand()rand(rb[1:3])
    @test A∧⋆(B, 𝑖) == (A⋅B)*𝑖
  end

  for i in 1:1000
    A = rand()rand(rb[1:5])
    B = rand()rand(rb[1:5])∧rand()rand(rb[1:5])
    @test A∧⋆(B) == (A⋅B)pseudoscalar(A)
    @test B∧⋆(A) == (B⋅A)pseudoscalar(A)
    C = rand()rand(rb[1:5])∧rand()rand(rb[1:5])∧rand()rand(rb[1:5])
    @test A∧⋆(C) == (A⋅C)pseudoscalar(A)
    @test C∧⋆(A) == (C⋅A)pseudoscalar(A)
    @test B∧⋆(C) == (B⋅C)pseudoscalar(A)
  end
end

