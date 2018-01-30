;;; Netlogo prototype implementation of the model
; Daily time step


;extensions [vid]
;type into observer
;vid:start-recorder
;vid:record-view ;; show the initial state
;repeat 30
;[ go
;  vid:record-view ]
;vid:save-recording "out.mp4"

globals
[
  wrong-C
  wrong-FONsum
  max_plantbiomass
]

breed [plants plant] ;#TODO: add others according to pollination syndrome later
breed [pollinators pollinator]

plants-own
[
  plant-biomass
  resource-FON-R ; use distance(xy) to calculate FON
  visited? ;used for reproduction
  juvenile?
  Fa
  C
  flowering?
]

pollinators-own
[
  poll-biomass
  dispersal-capability ;maximum distance of dispersal
  time-since-last-meal
  pollen ;for now, a float value that varies as the pollinator visits flowers
]

patches-own
[
  habitat?
  pollinator-resource
  plant-resource
  temperature
  FON-sum
  ;suitability
  n-FONs ;controls if FONs from different plants are added or the same FON just changes
  nonself-avgFa ; a way of avg other individuals FONs when there is overlap
]


to setup
  clear-all

  set-default-shape pollinators "my_bee"
  set-default-shape plants "flower"

  set wrong-C 0
  set max_plantbiomass 30

  ;;; Habitat: fragmented landscape ;;

   ;Set undisturbed landscape:
  ask patches
  [
    set habitat? true
    set pcolor 66
    set pollinator-resource 0
    set plant-resource 10
    set temperature random-normal mean-temperature 25
    set FON-sum 0
    set n-FONs 0
  ]

  ;;; Butterflies ;;;
  ;set them in teh habitat patches or
  create-pollinators n-pollinators
  [
    set size 3
    set color 2
    set poll-biomass random-normal 10 0.01
    set dispersal-capability maximum-dispersal
    set time-since-last-meal 999
    set pollen 0
  ]

  ;;; Plants ;;;
  create-plants n-plants
  [
    set size 5 ;size of the turtle shape, not the plant
    ;set color 14
    setxy random-pxcor random-pycor
    set plant-biomass 15
    set resource-FON-R (sqrt((plant-biomass ^ (2 / 3)) / pi) * 3) ;initial-FON-R
    set visited? true
    set juvenile? false ;they start mature

    ;setting up the FON: same as the plant-establish function
;    ask patches in-radius resource-FON-R
;      [
;        set n-FONs (n-FONs + 1) ;new individual, new FON
;        set FON-sum FON-sum + exp(-(distance myself))
;        ;set pcolor scale-color orange FON-sum 1 0
;        set pollinator-resource (pollinator-resource + 1)
;
;        ;procedure for ploting avg FON in the plot
;        if n-FONs >= 1 ;check if there is any FON in the patchs
;      [
;        ifelse n-FONs = 1
;        [set nonself-avgFa (FON-sum - exp(-(distance myself))) / n-FONs]
;        [set nonself-avgFa 0]
;      ]
;    ]

    ;setting simplified FON
    let FON exp(-(plant-biomass / max_plantbiomass))
    ask neighbors
    [
      let r (distance myself)
      set FON-sum (FON-sum +  FON * distance myself) ;adds
        ;decide on normalization rule for biomass: here, simple proportion
      set pcolor scale-color green FON-sum 1 0
    ]
  ]

  reset-ticks
end

to go

  ; Pollinators ;;;
  ask pollinators
  [
    search
    ;disperse ;TODO after a time without a meal, it will "search" farther
    pollinate
    poll-reproduce
    poll-die
  ]

  ;;; Plants ;;;
  ask plants
  [
    plant-establish ;recruitment is pop level. Individuals establish
    plant-grow
    flower
    plant-reproduce
    get-old
    plant-die
  ]

  ask patches ;just maintaining fragmented viasualization
  [
    ifelse habitat?
    [set pcolor 66]
    [set pcolor 36]
  ]

  ;;; DISTURBANCE every 1 years #TODO maybe even less
  if (ticks > 0 and remainder ticks 50 = 0 and ticks < 200) ;the disturbance stops, but the debt will still be paid
  [
    ask patches
    [set habitat? false]
     ;Set fragments
     ;TODO check how to set them all at once
    ask patches with [pxcor = -30 and pycor = 0]
    [
      ask patches in-radius (30 - (ticks / 100 + 1) * 3) ;as time goes on, the habitat fragments become smaller
      [
        set habitat? true
        set pcolor 66
        ;set pollinator-resource 0
        ;set plant-resource 10
        ;set temperature random-normal mean-temperature 5
        ;set FON-sum 0
        ;set n-FONs 0
      ]
    ]
    ask patches with [pxcor = 20 and pycor = -20]
    [
      ask patches in-radius (20 - (ticks / 100 + 1) * 3)
      [
        set habitat? true
        set pcolor 66
        ;set pollinator-resource 0
        ;set plant-resource 10
        ;set temperature random-normal mean-temperature 5
        ;set FON-sum 0
        ;set n-FONs 0
      ]
    ]
    ask patches with [pxcor = 25 and pycor = 25]
    [
      ask patches in-radius (25 - (ticks / 100 + 1) * 3)
      [
        set habitat? true
        set pcolor 66
;        set pollinator-resource 0
;        set plant-resource 10
;        set temperature random-normal mean-temperature 5
;        set FON-sum 0
;        set n-FONs 0
      ]
    ]

    ;Everything else is destroyed habitat
    ask patches with [habitat?  = false]
    [
      set temperature random-normal mean-temperature 25
      set n-FONs 0
      set pcolor 36
      set FON-sum 0
      set plant-resource  plant-resource - 5
    ]
  ]

  ask patches with [habitat? = true]
  [
    set plant-resource (plant-resource + (0.1 * (count neighbors with [habitat? = true])/(count neighbors))) ;a patch habitat benefits from being surrounded by others
  ]
  ; if I make this advantage dependent on the whole amount of habitat patches, it slows down a lot

  tick
  ;export-view (word "fonib" ticks ".png")
end

;;; POLLINATOR PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to search
  ;search nearest floral resoure
  let closestflower max-one-of neighbors [pollinator-resource] ;dont use uphill because bee stays put when there are no pollinator ressources around
  ifelse [pollinator-resource] of closestflower > pollinator-resource
  [
    move-to closestflower
    ;loose some biomass while searching
    set poll-biomass (poll-biomass - 0.3 ) ;#TODO: ajust to penalize dispersal
  ]
  [
    right (random 181) - 90
    let dist-search random 10
    forward dist-search
    ;loose some biomass while searching
    set poll-biomass (poll-biomass - 0.1 * dist-search) ;#TODO: ajust to penalize dispersal
  ]

end

to pollinate
  ; check for floral resources around
  ifelse poll-biomass > 10
  [
    set time-since-last-meal time-since-last-meal + 1
    stop
  ] ;too fat, no need
  [
    ifelse [pollinator-resource] of patch-here > 2 ;plant are dying too easy
      [
        ;Patch looses some floral resource
        ask patch-here
        [
          set pollinator-resource pollinator-resource - 1
        ]
        ;Plant looses biomass and "registers2 the visit
        ask (min-one-of plants [distance myself]) ;gets the closest plant
          [
            set plant-biomass (plant-biomass - 0.01) ; #TODO use MTE?
            set visited? true
          ]
        ;Pollinator gains biomass
        set poll-biomass (poll-biomass + 0.4) ; #TODO use MTE rates
        set time-since-last-meal 0
        ; and carries pollen ;TODO might not be necessary (this is supposed to be simple)
        set pollen (pollen + 0.1)
      ]
      [
      stop
      ]
    ]
end

to poll-reproduce
  if poll-biomass > 3 and (5 < random 11)
  [
    hatch 1 ;1 juvenile 'TODO make it more realistic: relevan pro Julia: eusocial and social might be different for the quantity of offspring AND eco-evol because socils are clones!
    [
      set poll-biomass random-normal 6 0.01
      set dispersal-capability maximum-dispersal
      set time-since-last-meal 999
      setxy random-xcor random-ycor
      if (not [habitat?] of patch-here)
      [die]
    ]
  ]
   ;loose some biomass
    set poll-biomass (poll-biomass - 0.2)
end

to poll-die
  if (ticks > 0 and random 11 > 5) or poll-biomass < 1 ;or (remainder [age] 14 0) #TODO densiy indeendent mortality
  [die]
end

;;; PLANT PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to plant-establish
  let suitability ([FON-sum] of patch-here < 0.8) ;suitable if total FON < 0.8

  if (juvenile? and suitability)
  [
    set juvenile? false

  ]
end

to plant-grow
  set plant-biomass plant-biomass + 0.05 * ([plant-resource] of patch-here) ;growth depends on local resources

  ask patches in-radius resource-FON-R
  [
    set plant-resource plant-resource - (0.001 * plant-resource)
  ]

  let FON exp(-(plant-biomass / max_plantbiomass))
  ask neighbors
  [
    let r (distance myself)
    set FON-sum (FON-sum +  FON * distance myself) ;adds
      ;decide on normalization rule for biomass: here, simple proportion
    set pcolor scale-color green FON-sum 1 0
  ]
end

to flower
  if remainder ticks 30 = 0
  [
    ask patches in-radius resource-FON-R
    [
      set pollinator-resource pollinator-resource + 2
    ]
  ]
end

to plant-reproduce
  if (plant-biomass > 10) and (ticks > 0 and remainder ticks 50 = 0) and (visited?)
  [
      hatch 1
      [
        ;plant-establish
        set size 5 ;size of the turtle shape, not the plant
        setxy random-xcor random-ycor ;#TODO: set a dispersal rate here
        set plant-biomass random-normal 10 3
        set resource-FON-R (sqrt((plant-biomass ^ (2 / 3)) / pi) * 3) ;initial-FON-R
        set visited? false
        set juvenile? true ;will establish a FON in next time step
    ]
   set plant-biomass plant-biomass - 1.5
    set visited? false ;will only reproduce again if visited
  ]
end

to get-old
  set plant-biomass (plant-biomass - 0.1)
end

to plant-die
  if plant-biomass < 5.5  or ([plant-resource] of patch-here <= 0) or 3 > random 1000
  [
    ask patches in-radius resource-FON-R
      [
        ifelse n-FONs > 0
        [
          set pollinator-resource (pollinator-resource - (pollinator-resource / n-FONs)) ;n-FONs in case it goes to 0
        ]
        [
          set pollinator-resource 0
        ]
          set FON-sum (FON-sum - exp(-(distance myself))) ; the total FON there should decrease
          set pcolor scale-color green FON-sum 10 0
          set n-FONs n-FONs - 1
      ]
    die
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
1070
871
-1
-1
12.0
1
10
1
1
1
0
0
0
1
-35
35
-35
35
1
1
1
ticks
15.0

BUTTON
33
773
106
806
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
113
774
176
807
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

INPUTBOX
54
147
163
207
maximum-dispersal
5.0
1
0
Number

INPUTBOX
55
80
162
140
mean-temperature
25.0
1
0
Number

INPUTBOX
107
10
191
70
n-plants
300.0
1
0
Number

INPUTBOX
21
10
105
70
n-pollinators
60.0
1
0
Number

MONITOR
1311
485
1469
530
NIL
count pollinators
17
1
11

PLOT
1209
232
1846
477
Populations monitoring
Time
Abundance
0.0
100.0
0.0
100.0
true
true
"" ""
PENS
"Pollinators" 1.0 0 -955883 true "" "plot count pollinators"
"Plants" 1.0 0 -13840069 true "" "plot count plants"

MONITOR
1209
484
1303
529
NIL
count plants
17
1
11

@#$#@#$#@
## WHAT IS IT?
(a general understanding of what the model is trying to show or explain)

This is a simplified version of the FONIB model that simulates the dynamics of interactions of metacommunities. Landscape resolution/specification is limited to the article juliano sent me.

## HOW IT WORKS

Plants will use the resources that are inside their zone of influence (ZOI). When zones of neighboring plants overlap in a patch, plants will compete for its resources with strength that is decreasingly proportional to the patch's distance to their stems. This results in their field-of-neighborhood (FON). Competition for pollinators is described in the same way: each plant has a ZOI for its floral resources and and will compete for eventual mutual pollinators. Pollinators are the equivalent to ressources here, but the funcionamento is not identical: i) the overlap for *pollination* FONs attracts more pollinators, which increases the chances of pollination for all species involved, whereas for *resources* FONs, it only decreases the availability of resources and ii) pollinators indirectly beneficial to reproduction (they are consuming plant biomass, but this is not a problem, because this is what such biomass is supposed to do - is this a relevant consideration?), but might also be detrimental, in case the visit is not effective.
In this first version of the model, abundances of an ecologically specialist and a generalist pollinator species and a strong and weak competitor plant species are monitored under different scenarios of temperature change

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

bee 2
true
0
Polygon -1184463 true false 195 150 105 150 90 165 90 225 105 270 135 300 165 300 195 270 210 225 210 165 195 150
Rectangle -16777216 true false 90 165 212 185
Polygon -16777216 true false 90 207 90 226 210 226 210 207
Polygon -16777216 true false 103 266 198 266 203 246 96 246
Polygon -6459832 true false 120 150 105 135 105 75 120 60 180 60 195 75 195 135 180 150
Polygon -6459832 true false 150 15 120 30 120 60 180 60 180 30
Circle -16777216 true false 105 30 30
Circle -16777216 true false 165 30 30
Polygon -7500403 true true 120 90 75 105 15 90 30 75 120 75
Polygon -16777216 false false 120 75 30 75 15 90 75 105 120 90
Polygon -7500403 true true 180 75 180 90 225 105 285 90 270 75
Polygon -16777216 false false 180 75 270 75 285 90 225 105 180 90
Polygon -7500403 true true 180 75 180 90 195 105 240 195 270 210 285 210 285 150 255 105
Polygon -16777216 false false 180 75 255 105 285 150 285 210 270 210 240 195 195 105 180 90
Polygon -7500403 true true 120 75 45 105 15 150 15 210 30 210 60 195 105 105 120 90
Polygon -16777216 false false 120 75 45 105 15 150 15 210 30 210 60 195 105 105 120 90
Polygon -16777216 true false 135 300 165 300 180 285 120 285

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
Polygon -8630108 true false 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -8630108 true false 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -955883 true false 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -955883 true false 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -7500403 true true 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -7500403 true true 135 90 30
Line -7500403 true 150 105 195 60
Line -7500403 true 150 105 105 60

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
Polygon -14835848 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
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

my_bee
true
0
Polygon -6459832 true false 120 75 105 90 105 120 135 150 165 150 195 120 195 90 180 75 120 75 120 75
Polygon -1184463 true false 152 149 105 165 75 195 67 211 74 234 85 252 100 264 116 276 134 286 150 285 167 285 182 278 206 260 220 242 226 218 225 195 195 165
Polygon -7500403 true true 151 69 119 74 105 90 90 75 90 60 90 45 103 33 120 45 135 30 150 30 165 30 180 45 197 32 210 45 210 60 210 75 195 90 180 75
Polygon -16777216 true false 70 185 74 171 223 172 224 186
Polygon -16777216 true false 67 211 71 226 224 226 225 211 67 211
Polygon -16777216 true false 91 257 106 269 195 269 211 255
Line -1 false 45 135 15 165
Line -1 false 15 165 30 195
Line -1 false 30 195 60 195
Line -1 false 90 165 135 90
Line -1 false 165 90 210 165
Line -1 false 210 165 240 195
Line -1 false 240 195 270 195
Line -1 false 285 165 255 135
Line -1 false 270 195 285 165
Line -1 false 255 135 165 90
Line -1 false 135 90 45 135
Line -1 false 60 195 90 165

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
NetLogo 6.0.2
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
