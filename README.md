# Quaternionic

A simple project to build a quaternion type in Julia.  I am unhappy with certain choices made in
[the most popular](https://github.com/JuliaGeometry/Quaternions.jl) quaternion package in Julia, and
I need a project to help me learn Julia itself, but I have a lot of experience programming
quaternions (including [this popular python package](https://github.com/moble/quaternion), and [this
newer and fancier package](https://github.com/moble/quaternionic)), so this seems like a good fit.
In particular, I am interested in understanding how a general geometric algebra would be coded, so I
am experimenting with static arrays, even though I can imagine that hard-coding the four components
would be faster.
