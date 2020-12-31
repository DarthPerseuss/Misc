globals
[
  tick-delta                           ;; how much we advance the tick counter this time through
  max-tick-delta
  left-particles right-particles       ;; particles in the left and right chambers
  avg-speed-cyan avg-energy-cyan       ;; left chamber averages
  avg-speed-magenta avg-energy-magenta ;; right chamber averages
  length-horizontal-surface-right           ;; the size of the wall surfaces that run horizontally - the top and bottom of the box
  length-horizontal-surface-left
  length-vertical-surface
  cyans magentas
  piston-kinetic-energy
  piston-vel
  piston-pos
]


breed [ particles particle ]
breed [ pistons piston ]
breed [ flashes flash ]

flashes-own [ birthday ]

pistons-own [ speed mass energy ]

particles-own
[
  speed mass energy          ;; particle info
  last-collision
  momentum-difference
  step-size
]


to setup
  clear-all
  set-default-shape particles "circle"
  set-default-shape flashes "square"
  set max-tick-delta 0.1073
  make-box
  make-piston
  set piston-vel 0
  set piston-pos init-piston-pos
  set length-vertical-surface  ( 2 * (box-height - 1) + 1)
  set length-horizontal-surface-right (box-width - init-piston-pos)
  set length-horizontal-surface-left  (init-piston-pos + box-width)
  make-particles
  update-variables
  reset-ticks
end


to update-variables
  set avg-speed-cyan      mean [speed] of cyans
  set avg-speed-magenta    mean [speed] of magentas
  set avg-energy-cyan     mean [energy] of cyans
  set avg-energy-magenta   mean [energy] of magentas

  ;; piston variables
  set piston-kinetic-energy (0.5 * piston-mass * (piston-vel ^ 2))
  set length-horizontal-surface-right (box-width - piston-pos)
  set length-horizontal-surface-left  (piston-pos + box-width)
end


to go
  if piston-pos < (- box-width + 2)
  [ user-message "The piston reached the left of the chamber. The simulation will stop."
    stop
  ]
  if piston-pos >= box-width - 2
  [ user-message "The piston reached the right of the chamber. The simulation will stop."
    stop
  ]

  ask particles [bounce]
  ask particles [ move ]
  ask particles [check-for-collision]
  move-piston
  ask particles [bounce]
  tick-advance tick-delta
  if floor ticks > floor (ticks - tick-delta)
  [
     update-variables
     update-plots
  ]
  calculate-tick-delta

  ask flashes with [ticks - birthday > 0.4 or pcolor = black]
    [ die ]
  display
end


to bounce  ;; particle procedure
  ;; get the coordinates of the patch we'll be on if we go forward 1
  let new-patch patch-ahead 1
  let new-px [pxcor] of new-patch
  let new-py [pycor] of new-patch
  ;; if we're not about to hit a wall, we don't need to do any further checks

  if ([pcolor] of new-patch != yellow and [pcolor] of new-patch != orange)
    [ stop ]
  ;; if hitting left or right wall, reflect heading around x axis
  if (abs new-px = box-width)
    [ set heading (- heading)
  ;;  if the particle is hitting a vertical wall, only the horizontal component of the speed
  ;;  vector can change.  The change in velocity for this component is 2 * the speed of the particle,
  ;;  due to the reversing of direction of travel from the collision with the wall
  ;;  if hitting top or bottom wall, reflect heading around y axis
  ]
  ;;  if the particle is hitting a horizontal wall, only the vertical component of the speed
  ;;  vector can change.  The change in velocity for this component is 2 * the speed of the particle,
  ;;  due to the reversing of direction of travel from the collision with the wall
  ;;  every time a particle hits the wall, it produces a short-living "flash" so assist in visualization
  if (abs new-py = box-height)
  [
   set heading (180 - heading)
  ]
  ;; hitting piston

  if ((new-px = [pxcor] of one-of pistons and ((speed * dx) > piston-vel or (speed * dx) < (- piston-vel))) or [pcolor] of new-patch = orange)
  [
    ;
    ;set momentum-difference momentum-difference + ((dx * 2 * mass * speed) / length-vertical-surface)   ;; no abs here, cz I think the direction is important
    exchange-energy-with-piston
  ]


  ask particles [safety-check]

  ask patch new-px new-py
    [
      sprout-flashes 1 [
        set color pcolor - 2
        set birthday ticks
        set heading 0
      ]
  ]
end


to safety-check
  if shade-of? color cyan and pcolor = orange
  [
    if heading <= 180
    [ set heading (heading + 360 - 2 * heading) ]
  ]
  if shade-of? color magenta and pcolor = orange
  [
    if heading > 180
      [ set heading (heading + 360 - 2 * heading)]
  ]

end

to move
    if patch-ahead (speed * tick-delta) != patch-here
      [ set last-collision nobody ]
    jump (speed * tick-delta)
end


to calculate-tick-delta
  ;; tick-delta is calculated in such way that even the fastest
  ;; particle will jump at most 1 patch length in a tick. As
  ;; particles jump (speed * tick-delta) at every tick, making
  ;; tick length the inverse of the speed of the fastest particle
  ;; (1/max speed) assures that. Having each particle advance at most
  ;; one patch-length is necessary for it not to "jump over" a wall.
    ifelse any? particles with [speed > 0]
    [ set tick-delta min list
      (1 / (ceiling max (sentence ([speed] of particles) ([speed] of one-of pistons))))
      max-tick-delta ]
    [ set tick-delta max-tick-delta ]
end


to check-for-collision
  if count other particles-here = 1
  [
    ;; the following conditions are imposed on collision candidates:
    ;;   1. they must have a lower who number than my own, because collision
    ;;      code is asymmetrical: it must always happen from the point of view
    ;;      of just one particle.
    ;;   2. they must not be the same particle that we last collided with on
    ;;      this patch, so that we have a chance to leave the patch after we've
    ;;      collided with someone.
    let candidate one-of other particles-here with
      [who < [who] of myself and myself != last-collision]
    ;; we also only collide if one of us has non-zero speed. It's useless
    ;; (and incorrect, actually) for two particles with zero speed to collide.
    if (candidate != nobody) and (speed > 0 or [speed] of candidate > 0)
    [
      collide-with candidate
      set last-collision candidate
      ask candidate [ set last-collision myself ]
    ]
  ]
end


to collide-with [ other-particle ] ;; particle procedure
  ;;; PHASE 1: initial setup

  ;; for convenience, grab some quantities from other-particle
  let mass2 [mass] of other-particle
  let speed2 [speed] of other-particle
  let heading2 [heading] of other-particle

  ;; since particles are modeled as zero-size points, theta isn't meaningfully
  ;; defined. we can assign it randomly without affecting the model's outcome.
  let theta (random-float 360)



  ;;; PHASE 2: convert velocities to theta-based vector representation

  ;; now convert my velocity from speed/heading representation to components
  ;; along theta and perpendicular to theta
  let v1t (speed * cos (theta - heading))
  let v1l (speed * sin (theta - heading))

  ;; do the same for other-particle
  let v2t (speed2 * cos (theta - heading2))
  let v2l (speed2 * sin (theta - heading2))



  ;;; PHASE 3: manipulate vectors to implement collision

  ;; compute the velocity of the system's center of mass along theta
  let vcm (((mass * v1t) + (mass2 * v2t)) / (mass + mass2) )

  ;; now compute the new velocity for each particle along direction theta.
  ;; velocity perpendicular to theta is unaffected by a collision along theta,
  ;; so the next two lines actually implement the collision itself, in the
  ;; sense that the effects of the collision are exactly the following changes
  ;; in particle velocity.
  set v1t (2 * vcm - v1t)
  set v2t (2 * vcm - v2t)



  ;;; PHASE 4: convert back to normal speed/heading

  ;; now convert my velocity vector into my new speed and heading
  set speed sqrt ((v1t ^ 2) + (v1l ^ 2))
  set energy (0.5 * mass * speed * speed)
  ;; if the magnitude of the velocity vector is 0, atan is undefined. but
  ;; speed will be 0, so heading is irrelevant anyway. therefore, in that
  ;; case we'll just leave it unmodified.
  if v1l != 0 or v1t != 0
    [ set heading (theta - (atan v1l v1t)) ]

  ;; and do the same for other-particle
  ask other-particle [
    set speed sqrt ((v2t ^ 2) + (v2l ^ 2))
    set energy (0.5 * mass * (speed ^ 2))
    if v2l != 0 or v2t != 0
      [ set heading (theta - (atan v2l v2t)) ]
  ]

    ;; PHASE 5: final updates

  ;; now recolor, since color is based on quantities that may have changed

  recolor
  ask other-particle
    [ recolor ]
end

to recolor  ;; particle procedure
  let values [ speed ] of breed
  let lower-limit mean values - 3 * standard-deviation values
  let upper-limit mean values + 3 * standard-deviation values
  set color scale-color color speed (lower-limit - 1) (upper-limit + 1)
end



to calculate-pressure-left



end


to calculate-pressure-right




end


to make-box
  ask patches with [((abs pxcor = box-width) and  (abs pycor <= box-height)) or
                    ((abs pycor = box-height) and (abs pxcor <= box-width))]
    [ set pcolor yellow ]
end


;; creates initial particles
to make-particles
  create-particles num-cyans [
    set speed cyan-init-speed
    set mass cyan-mass
    random-position-left
    set color cyan
  ]
  create-particles num-magentas [
    set speed magenta-init-speed
    set mass magenta-mass
    random-position-right
    set shape "circle"
    set color magenta
  ]
  set cyans particles with [color = cyan]
  set magentas particles with [color = magenta]
  ask particles
  [
    set energy (0.5 * mass * speed * speed)
    ;; make their graphical size equal to the cube root of their mass
    set size  mass ^ .1
    set last-collision nobody
  ]
  calculate-tick-delta
end

;; place particle at random location inside the box.
to random-position-left  ;; particle procedure
  setxy ( - box-width + 1 + random-float (box-width - init-piston-pos - 1.5) )
        ( - box-height + 1 + 2 * random-float (box-height - 1.5) )
end

to random-position-right
    setxy ( init-piston-pos + 1 + random-float (box-width - init-piston-pos - 1.5))
          (- box-height + 1 + 2 * random-float (box-height - 1.5))
end
;; ------ Piston ----------
to make-piston
  ask patches with [pycor >= (- box-height + 1) and pycor <= box-height - 1 and pxcor = init-piston-pos]
  [ sprout-pistons 1
    [ set color orange
      set pcolor color
      hide-turtle
    ]
  ]
end

to move-piston
  let old-piston-vel piston-vel
  let movement-amount (old-piston-vel * tick-delta)
  ;; Setting the pcolor makes the piston look like a wall to the particles.
  ask pistons
  [ set pcolor black
    while [(piston-vel * tick-delta) >= 1.0]
      [ calculate-tick-delta ]
    if piston-vel > 0 [
      ifelse piston-pos + movement-amount <=  (box-width - 1)
      [
        set heading 90
        fd movement-amount
        set piston-pos ([xcor] of one-of pistons)
      ]
      [
        set xcor (box-width - 1)
        set piston-pos (box-width - 1)
        if piston-vel > 0
            [ set piston-vel 0 ]
      ]
    ]

    if piston-vel < 0 [
      ifelse piston-pos + movement-amount >=  (- box-width + 1)
      [
        set heading 270
        bk movement-amount
        set piston-pos ([xcor] of one-of pistons)
      ]
      [
        set xcor (- box-width + 1)
        set piston-pos (- box-width + 1)
        if piston-vel < 0
            [ set piston-vel 0 ]
      ]
    ]

    set speed piston-vel ;; just used for tick-delta calculations
    set pcolor color

    if (piston-vel < 0)  ;; when piston going to the left
      [ if (any? particles-here with [(speed * dx) > piston-vel])
        [ ;; only bounce particles that are moving down slower than the piston
          ;; faster ones should outrun it
          ask particles-here with [(speed * dx) > piston-vel]
          [
            ;;  if the particle is hitting the piston, only the vertical component of the speed
            ;;  vector can change.  The change in velocity for this component is 2 * the speed of the particle,
            ;; due to the reversing of direction of travel from the collision with the wall
            ;; make sure that each particle finishes exchanging energy before any others can
            exchange-energy-with-piston
          ]
        ]
    ]

    if (piston-vel > 0)  ;; when piston going to the right
      [ if (any? particles-here with [(speed * dx) < piston-vel])
        [; piston moving to the right
          ask particles-here with [(speed * dx) < piston-vel]
          [
            exchange-energy-with-piston
          ]
        ]
    ]
  ]
end


to exchange-energy-with-piston  ;; particle procedure -- piston and particle exchange energy
  let vx (speed * dx)         ;;only along x-axis
  let vy (speed * dy)         ;;only along y-axis
  let old-vx vx
  let old-piston-vel piston-vel
  set piston-vel ((((piston-mass - mass) / (piston-mass + mass)) * old-piston-vel) +
                  (((2 * mass) / (piston-mass + mass)) * old-vx))
  set vx ((((2 * piston-mass) / (piston-mass + mass)) * old-piston-vel) -
         (((piston-mass - mass) / (piston-mass + mass)) * old-vx))
  set speed (sqrt ((vx ^ 2) + (vy ^ 2)))
  set energy (0.5 * mass * (speed  ^ 2))
  set heading atan vx vy  ;; horizontal collision
  recolor
end


@#$#@#$#@
GRAPHICS-WINDOW
501
18
1055
573
-1
-1
6.0
1
10
1
1
1
0
0
0
1
-45
45
-45
45
1
1
1
ticks
30.0

BUTTON
27
17
90
50
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
28
56
91
89
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
3
277
175
310
num-cyans
num-cyans
0
500
185.0
1
1
NIL
HORIZONTAL

SLIDER
3
318
175
351
cyan-init-speed
cyan-init-speed
0
100
49.9
0.1
1
NIL
HORIZONTAL

SLIDER
4
361
176
394
cyan-mass
cyan-mass
0
50
41.7
.1
1
NIL
HORIZONTAL

SLIDER
225
278
397
311
num-magentas
num-magentas
0
500
185.0
1
1
NIL
HORIZONTAL

SLIDER
226
320
398
353
magenta-init-speed
magenta-init-speed
0
100
50.0
.1
1
NIL
HORIZONTAL

SLIDER
226
362
398
395
magenta-mass
magenta-mass
0
100
38.2
.1
1
NIL
HORIZONTAL

SLIDER
225
91
397
124
init-piston-pos
init-piston-pos
-45
45
0.0
1
1
NIL
HORIZONTAL

SLIDER
225
17
397
50
box-width
box-width
1
45
43.0
1
1
NIL
HORIZONTAL

SLIDER
225
54
397
87
box-height
box-height
1
45
32.0
1
1
NIL
HORIZONTAL

SLIDER
225
129
397
162
piston-mass
piston-mass
20
1000
1000.0
20
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
