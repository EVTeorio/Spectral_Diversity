Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(40, 45, 52)
$muted = [System.Drawing.Color]::FromArgb(105, 112, 121)
$grid = [System.Drawing.Color]::FromArgb(215, 220, 226)
$green = [System.Drawing.Color]::FromArgb(38, 122, 80)
$blue = [System.Drawing.Color]::FromArgb(49, 104, 176)
$amber = [System.Drawing.Color]::FromArgb(188, 119, 34)
$red = [System.Drawing.Color]::FromArgb(184, 75, 68)
$teal = [System.Drawing.Color]::FromArgb(34, 137, 145)
$purple = [System.Drawing.Color]::FromArgb(114, 90, 168)

function Font($size, $style = [System.Drawing.FontStyle]::Regular) {
  return [System.Drawing.Font]::new($fontFamily, [single]$size, $style)
}

function Brush($color) { return [System.Drawing.SolidBrush]::new($color) }
function PenC($color, $width = 2) {
  if ($color -is [array]) { $color = $color[0] }
  if ($width -is [array]) { $width = $width[0] }
  return [System.Drawing.Pen]::new($color, [single]$width)
}

function Save-Figure($name, [scriptblock]$draw, $width = 2400, $height = 1600) {
  $bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  $g = [System.Drawing.Graphics]::FromImage($bmp)
  $g.Clear([System.Drawing.Color]::Transparent)
  $g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
  $g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit
  & $draw $g $width $height
  $path = Join-Path $outDir $name
  $bmp.Save($path, [System.Drawing.Imaging.ImageFormat]::Png)
  $g.Dispose()
  $bmp.Dispose()
  return $path
}

function Draw-Text($g, $text, $x, $y, $size = 30, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $f = Font $size $style
  $b = Brush $color
  $g.DrawString($text, $f, $b, [single]$x, [single]$y)
  $f.Dispose(); $b.Dispose()
}

function Draw-CenteredText($g, $text, $x, $y, $w, $h, $size = 28, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $f = Font $size $style
  $b = Brush $color
  $sf = [System.Drawing.StringFormat]::new()
  $sf.Alignment = [System.Drawing.StringAlignment]::Center
  $sf.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $f, $b, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $sf)
  $sf.Dispose(); $f.Dispose(); $b.Dispose()
}

function Draw-LabelBox($g, $x, $y, $w, $h, $title, $body, $color) {
  $fill = [System.Drawing.Color]::FromArgb(235, $color.R, $color.G, $color.B)
  $path = [System.Drawing.Drawing2D.GraphicsPath]::new()
  $r = 16
  $path.AddArc($x, $y, $r, $r, 180, 90)
  $path.AddArc($x + $w - $r, $y, $r, $r, 270, 90)
  $path.AddArc($x + $w - $r, $y + $h - $r, $r, $r, 0, 90)
  $path.AddArc($x, $y + $h - $r, $r, $r, 90, 90)
  $path.CloseFigure()
  $g.FillPath((Brush $fill), $path)
  $g.DrawPath((PenC $color 3), $path)
  Draw-Text $g $title ($x + 28) ($y + 22) 31 ([System.Drawing.Color]::White) ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g $body ($x + 28) ($y + 78) 25 ([System.Drawing.Color]::White) ([System.Drawing.FontStyle]::Bold)
  $path.Dispose()
}

function Draw-Axes($g, $left, $top, $right, $bottom, $title, $xlabel, $ylabel, $ymin, $ymax) {
  $lo = [double]$ymin
  $hi = [double]$ymax
  Draw-Text $g $title $left ($top - 82) 36 $ink ([System.Drawing.FontStyle]::Bold)
  $axisPen = PenC $ink 3
  $gridPen = PenC $grid 1
  $g.DrawLine($axisPen, $left, $bottom, $right, $bottom)
  $g.DrawLine($axisPen, $left, $top, $left, $bottom)
  for ($i = 0; $i -le 5; $i++) {
    $v = $lo + (($hi - $lo) * $i / 5)
    $y = $bottom - (($v - $lo) / ($hi - $lo)) * ($bottom - $top)
    $g.DrawLine($gridPen, $left, $y, $right, $y)
    Draw-Text $g ("{0:0.00}" -f $v) ($left - 88) ($y - 18) 22 $muted
  }
  Draw-CenteredText $g $xlabel (($left + $right) / 2 - 240) ($bottom + 80) 480 42 26 $ink
  $state = $g.Save()
  $g.TranslateTransform($left - 150, (($top + $bottom) / 2 + 180))
  $g.RotateTransform(-90)
  Draw-CenteredText $g $ylabel 0 0 360 42 26 $ink
  $g.Restore($state)
  $axisPen.Dispose(); $gridPen.Dispose()
}

function Scale-Y($value, $top, $bottom, $ymin, $ymax) {
  $v = [double]$value
  $lo = [double]$ymin
  $hi = [double]$ymax
  return $bottom - (($v - $lo) / ($hi - $lo)) * ($bottom - $top)
}

function Draw-LineSeries($g, $left, $top, $right, $bottom, $scales, $series, $ymin, $ymax) {
  $plotW = $right - $left
  $xPos = @{}
  for ($i = 0; $i -lt $scales.Count; $i++) {
    $xPos[$scales[$i]] = $left + ($plotW * ($i + 0.5) / $scales.Count)
    Draw-Text $g $scales[$i].Replace("m", " m") ($xPos[$scales[$i]] - 28) ($bottom + 22) 24 $ink
  }
  foreach ($s in $series) {
    $pen = PenC ($s.Color) 5
    $brush = Brush ($s.Color)
    $prev = $null
    foreach ($scale in $scales) {
      if (-not $s.Values.ContainsKey($scale)) { continue }
      $x = $xPos[$scale]
      $y = Scale-Y $s.Values[$scale] $top $bottom $ymin $ymax
      if ($null -ne $prev) { $g.DrawLine($pen, $prev[0], $prev[1], $x, $y) }
      $g.FillEllipse($brush, $x - 10, $y - 10, 20, 20)
      Draw-Text $g ("{0:0.000}" -f [double]$s.Values[$scale]) ($x - 36) ($y - 46) 21 $s.Color
      $prev = @($x, $y)
    }
    $pen.Dispose(); $brush.Dispose()
  }
}

function Draw-Legend($g, $items, $x, $y) {
  $yy = $y
  foreach ($item in $items) {
    $pen = PenC ($item.Color) 5
    $g.DrawLine($pen, $x, $yy + 16, $x + 52, $yy + 16)
    $g.FillEllipse((Brush ($item.Color)), $x + 18, $yy + 6, 20, 20)
    Draw-Text $g $item.Label ($x + 72) $yy 24 $ink
    $pen.Dispose()
    if ($item.Label -match "`n") { $yy += 82 } else { $yy += 46 }
  }
}

function Label-Spectral($metric) {
  switch ($metric) {
    "spec_spca_alpha" { "standardized PCA alpha hull" }
    "spec_spca_mean" { "standardized PCA mean distance" }
    "spec_spca_rao" { "standardized PCA Rao's Q" }
    "spec_alpha" { "raw PCA alpha hull" }
    "spec_pca_mean" { "raw PCA mean distance" }
    "spec_rao_q" { "raw PCA Rao's Q" }
    default { $metric }
  }
}

function Label-Bio($metric) {
  switch ($metric) {
    "phy_faith" { "Faith's PD" }
    "phy_rao" { "phylogenetic Rao's Q" }
    "phy_afaith" { "abundance-weighted Faith's PD" }
    "sp_rich" { "species richness" }
    "sp_shannon" { "Shannon diversity" }
    "sp_simpson" { "Simpson diversity" }
    "sp_even" { "species evenness" }
    default { $metric }
  }
}

function Fmt-R($x) { return ("{0:0.000}" -f [double]$x) }
function Fmt-R2($x) { return ("{0:0.000}" -f [double]$x) }
function Fmt-P($p) {
  $d = [double]$p
  if ($d -eq 0) { return "<0.001" }
  if ($d -lt 0.001) { return "<0.001" }
  return ("{0:0.003}" -f $d)
}

$all = Import-Csv (Join-Path $root "reports/tables/final_research_direction/all_spectral_biodiversity_relationships.csv")
$spectralPairs = Import-Csv (Join-Path $root "reports/tables/final_research_direction/spectral_metric_pairwise_relationships.csv")
$bioPairs = Import-Csv (Join-Path $root "reports/tables/final_research_direction/biodiversity_metric_pairwise_relationships.csv")
$moran = Import-Csv (Join-Path $root "reports/tables/final_research_direction/spatial_moran_diagnostics.csv")
$trans = Import-Csv (Join-Path $root "reports/tables/final_research_direction/transformation_best_relationship_summary.csv")
$rawPc = Import-Csv (Join-Path $root "reports/tables/pc1_mean_reflectance_correlation/50m_regular_PCA_pc_mean_reflectance_correlation_summary.csv")
$stdPc = Import-Csv (Join-Path $root "reports/tables/pc1_mean_reflectance_correlation/50m_standardized_PCA_pc_mean_reflectance_correlation_summary.csv")

$scales = @("10m", "20m", "50m")

Save-Figure "01_study_site_quadrat_scales_workflow.png" {
  param($g, $w, $h)
  Draw-Text $g "Study site, quadrat scales, and workflow" 100 70 44 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g "Paint Rock Forest Dynamics Plot | drone-acquired VNIR hyperspectral imagery | canopy-weighted biodiversity" 100 130 26 $muted
  Draw-LabelBox $g 120 240 430 190 "20 ha forest plot" "Mixed hardwood-pine`nNorthern Alabama" $green
  Draw-LabelBox $g 640 240 470 190 "2024 forest census" "4,007 canopy trees`n51 species" $blue
  Draw-LabelBox $g 1220 240 450 190 "VNIR imagery" "June 2024 drone data`nshadow-screened spectra" $teal
  Draw-LabelBox $g 1780 240 500 190 "PCA spectral space" "Raw PCA`nvector-normalized PCA" $purple
  $arrow = PenC $muted 5
  $arrow.EndCap = [System.Drawing.Drawing2D.LineCap]::ArrowAnchor
  $g.DrawLine($arrow, 570, 335, 620, 335)
  $g.DrawLine($arrow, 1130, 335, 1200, 335)
  $g.DrawLine($arrow, 1690, 335, 1760, 335)
  $arrow.Dispose()

  Draw-Text $g "Shared quadrat grains" 145 570 34 $ink ([System.Drawing.FontStyle]::Bold)
  $boxPen = PenC $ink 4
  $colors = @($green, $blue, $purple)
  $labels = @("10 m", "20 m", "50 m")
  for ($i = 0; $i -lt 3; $i++) {
    $x = 150 + $i * 235
    $s = 145 + $i * 45
    $g.DrawRectangle((PenC ($colors[$i]) 5), $x, 650, $s, $s)
    Draw-CenteredText $g $labels[$i] $x (650 + $s + 18) $s 44 26 $colors[$i] ([System.Drawing.FontStyle]::Bold)
  }
  $boxPen.Dispose()

  Draw-LabelBox $g 790 610 620 230 "Biodiversity" "Species richness, Shannon diversity,`nSimpson diversity, species evenness,`nFaith's PD, phylogenetic Rao's Q,`nabundance-weighted Faith's PD" $blue
  Draw-LabelBox $g 1590 610 620 230 "Spectral metrics" "standardized PCA alpha hull,`nstandardized PCA mean distance,`nstandardized PCA Rao's Q" $green
  $arrow = PenC $muted 5
  $arrow.EndCap = [System.Drawing.Drawing2D.LineCap]::ArrowAnchor
  $g.DrawLine($arrow, 1430, 755, 1570, 755)
  $arrow.Dispose()

  Draw-LabelBox $g 500 1040 1400 230 "Analysis focus" "Scale, metric choice, and preprocessing`nshape spectral-biodiversity relationships." $amber
} | Out-Null

$spSeries = @()
$spWanted = @(
  @{ X = "spec_spca_alpha"; Y = "spec_spca_mean"; Label = "alpha hull vs. mean distance"; Color = $green },
  @{ X = "spec_spca_mean"; Y = "spec_spca_rao"; Label = "mean distance vs. Rao's Q"; Color = $blue },
  @{ X = "spec_spca_alpha"; Y = "spec_spca_rao"; Label = "alpha hull vs. Rao's Q"; Color = $amber }
)
foreach ($want in $spWanted) {
  $vals = @{}
  foreach ($scale in $scales) {
    $row = $spectralPairs | Where-Object { $_.scale -eq $scale -and (($_.x_metric -eq $want.X -and $_.y_metric -eq $want.Y) -or ($_.x_metric -eq $want.Y -and $_.y_metric -eq $want.X)) } | Select-Object -First 1
    if ($row) { $vals[$scale] = [double]$row.pearson_r }
  }
  $spSeries += [pscustomobject]@{ Label = $want.Label; Color = $want.Color; Values = $vals }
}
Save-Figure "02_standardized_pca_spectral_metric_relationships_by_scale.png" {
  param($g, $w, $h)
  $left = 260; $top = 210; $right = 1760; $bottom = 1240
  Draw-Axes $g $left $top $right $bottom "Spectral metric relationships by scale" "Quadrat grain" "Pearson r" 0 1
  Draw-LineSeries $g $left $top $right $bottom $scales $spSeries 0 1
  Draw-Legend $g $spSeries 1840 300
  Draw-Text $g "Metrics use standardized PCA space." 1840 520 23 $muted
} | Out-Null

$bioSeries = @()
$bioWanted = @(
  @{ X = "sp_shannon"; Y = "sp_simpson"; Label = "Shannon vs. Simpson"; Color = $blue },
  @{ X = "sp_rich"; Y = "sp_shannon"; Label = "richness vs. Shannon"; Color = $green },
  @{ X = "phy_faith"; Y = "sp_rich"; Label = "Faith's PD vs. species richness"; Color = $amber },
  @{ X = "phy_rao"; Y = "phy_afaith"; Label = "phylogenetic Rao's Q vs.`nabundance-weighted Faith's PD"; Color = $purple },
  @{ X = "sp_shannon"; Y = "phy_afaith"; Label = "Shannon diversity vs.`nabundance-weighted Faith's PD"; Color = $red }
)
foreach ($want in $bioWanted) {
  $vals = @{}
  foreach ($scale in $scales) {
    $row = $bioPairs | Where-Object { $_.scale -eq $scale -and (($_.x_metric -eq $want.X -and $_.y_metric -eq $want.Y) -or ($_.x_metric -eq $want.Y -and $_.y_metric -eq $want.X)) } | Select-Object -First 1
    if ($row) { $vals[$scale] = [double]$row.pearson_r }
  }
  $bioSeries += [pscustomobject]@{ Label = $want.Label; Color = $want.Color; Values = $vals }
}
Save-Figure "03_biodiversity_metric_relationships_by_scale.png" {
  param($g, $w, $h)
  $left = 260; $top = 210; $right = 1700; $bottom = 1240
  Draw-Axes $g $left $top $right $bottom "Biodiversity metric relationships by scale" "Quadrat grain" "Pearson r" -0.2 1
  Draw-LineSeries $g $left $top $right $bottom $scales $bioSeries -0.2 1
  Draw-Legend $g $bioSeries 1780 260
} | Out-Null

$bioMetrics = @("phy_afaith", "phy_rao", "phy_faith", "sp_rich", "sp_shannon", "sp_simpson", "sp_even")
$specMetrics = @("spec_spca_alpha", "spec_spca_mean", "spec_spca_rao")
Save-Figure "04_standardized_pca_spectral_biodiversity_correlations_by_scale.png" {
  param($g, $w, $h)
  Draw-Text $g "Spectral-biodiversity correlations by scale" 100 70 44 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g "Pearson r for standardized PCA spectral heterogeneity metrics" 100 130 27 $muted
  $cellW = 145; $cellH = 78
  $startX = 590; $startY = 250
  for ($p = 0; $p -lt $scales.Count; $p++) {
    $x0 = $startX + $p * 610
    Draw-Text $g $scales[$p].Replace("m"," m") ($x0 + 170) 200 32 $ink ([System.Drawing.FontStyle]::Bold)
    for ($c = 0; $c -lt $specMetrics.Count; $c++) {
      Draw-CenteredText $g (Label-Spectral $specMetrics[$c]).Replace("standardized PCA ","").Replace("mean distance","mean`ndistance").Replace("alpha hull","alpha`nhull") ($x0 + $c * $cellW) 250 $cellW 90 21 $ink
    }
    for ($r = 0; $r -lt $bioMetrics.Count; $r++) {
      if ($p -eq 0) { Draw-Text $g (Label-Bio $bioMetrics[$r]) 80 (360 + $r * $cellH + 18) 22 $ink }
      for ($c = 0; $c -lt $specMetrics.Count; $c++) {
        $row = $all | Where-Object { $_.scale -eq $scales[$p] -and $_.x_metric -eq $bioMetrics[$r] -and $_.y_metric -eq $specMetrics[$c] } | Select-Object -First 1
        $rv = if ($row) { [double]$row.pearson_r } else { 0.0 }
        if ($rv -ge 0) {
          $intensity = [Math]::Min(235, [Math]::Max(25, [int](35 + 380 * [Math]::Abs($rv))))
          $col = [System.Drawing.Color]::FromArgb(45 + $intensity / 3, 230 - $intensity / 3, 245 - $intensity / 2)
        } else {
          $intensity = [Math]::Min(220, [Math]::Max(20, [int](60 + 380 * [Math]::Abs($rv))))
          $col = [System.Drawing.Color]::FromArgb(250, 240 - $intensity / 2, 230 - $intensity / 2)
        }
        $rect = [System.Drawing.Rectangle]::new($x0 + $c * $cellW, 360 + $r * $cellH, $cellW - 8, $cellH - 8)
        $g.FillRectangle((Brush $col), $rect)
        $g.DrawRectangle((PenC ([System.Drawing.Color]::FromArgb(245,245,245)) 2), $rect)
        Draw-CenteredText $g (Fmt-R $rv) $rect.X $rect.Y $rect.Width $rect.Height 25 $ink
      }
    }
  }
} | Out-Null

$scaleSeries = @()
foreach ($want in @(
  @{ Y = "spec_spca_alpha"; Label = "standardized PCA`nalpha hull"; Color = $green },
  @{ Y = "spec_spca_mean"; Label = "standardized PCA`nmean distance"; Color = $blue },
  @{ Y = "spec_spca_rao"; Label = "standardized PCA`nRao's Q"; Color = $amber }
)) {
  $vals = @{}
  foreach ($scale in $scales) {
    $row = $all | Where-Object { $_.scale -eq $scale -and $_.x_metric -eq "phy_afaith" -and $_.y_metric -eq $want.Y } | Select-Object -First 1
    if ($row) { $vals[$scale] = [double]$row.pearson_r }
  }
  $scaleSeries += [pscustomobject]@{ Label = $want.Label; Color = $want.Color; Values = $vals }
}
Save-Figure "05_abundance_weighted_faith_pd_scale_response.png" {
  param($g, $w, $h)
  $left = 260; $top = 210; $right = 1720; $bottom = 1240
  Draw-Axes $g $left $top $right $bottom "Scale response for abundance-weighted Faith's PD" "Quadrat grain" "Pearson r" 0 0.55
  Draw-LineSeries $g $left $top $right $bottom $scales $scaleSeries 0 0.55
  Draw-Legend $g $scaleSeries 1800 320
} | Out-Null

$pcRows = @()
foreach ($row in $rawPc + $stdPc) {
  if ($row.pc_axis -eq "PC1") {
    $pcRows += $row
  }
}
Save-Figure "06_raw_vs_vector_normalized_pca_brightness_relationship.png" {
  param($g, $w, $h)
  Draw-Text $g "Raw PCA versus vector-normalized PCA brightness relationship" 100 70 42 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g "50 m analysis; PC1 correlated with mean reflectance" 100 130 27 $muted
  $left = 260; $top = 250; $right = 1780; $bottom = 1230
  Draw-Axes $g $left $top $right $bottom "" "PCA basis and analysis level" "R2" 0 0.9
  $bars = @(
    @{ Label = "raw PCA`npixel"; Value = ($pcRows | Where-Object {$_.pca_basis -eq "regular_PCA" -and $_.analysis_level -eq "pixel_level"} | Select-Object -First 1).r_squared; Color = $amber },
    @{ Label = "raw PCA`nquadrat mean"; Value = ($pcRows | Where-Object {$_.pca_basis -eq "regular_PCA" -and $_.analysis_level -eq "quadrat_mean_level"} | Select-Object -First 1).r_squared; Color = $amber },
    @{ Label = "vector-normalized PCA`npixel"; Value = ($pcRows | Where-Object {$_.pca_basis -eq "standardized_PCA" -and $_.analysis_level -eq "pixel_level"} | Select-Object -First 1).r_squared; Color = $green },
    @{ Label = "vector-normalized PCA`nquadrat mean"; Value = ($pcRows | Where-Object {$_.pca_basis -eq "standardized_PCA" -and $_.analysis_level -eq "quadrat_mean_level"} | Select-Object -First 1).r_squared; Color = $green }
  )
  for ($i = 0; $i -lt $bars.Count; $i++) {
    $x = $left + 120 + $i * 330
    $v = [double]$bars[$i].Value
    $y = Scale-Y $v $top $bottom 0 0.9
    $g.FillRectangle((Brush ($bars[$i].Color)), $x, $y, 170, $bottom - $y)
    Draw-CenteredText $g ("{0:0.000}" -f $v) ($x - 20) ($y - 55) 210 44 25 $bars[$i].Color ([System.Drawing.FontStyle]::Bold)
    Draw-CenteredText $g $bars[$i].Label ($x - 60) ($bottom + 30) 290 80 22 $ink
  }
} | Out-Null

$moranKeep = $moran | Where-Object {
  ($_.diagnostic_type -eq "variable" -and $_.variable -in @("spec_spca_alpha","spec_spca_mean","spec_spca_rao","phy_afaith","phy_rao","sp_shannon")) -or
  ($_.diagnostic_type -eq "residual" -and $_.model -in @("spec_spca_alpha ~ phy_afaith","spec_spca_alpha ~ phy_rao","spec_spca_alpha ~ sp_shannon","spec_spca_mean ~ phy_afaith","spec_spca_mean ~ phy_rao","spec_spca_mean ~ sp_shannon"))
}
Save-Figure "07_spatial_moran_i_diagnostics.png" {
  param($g, $w, $h)
  Draw-Text $g "Spatial autocorrelation diagnostics" 100 70 44 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g "Moran's I for primary variables and residual relationships" 100 130 27 $muted
  $left = 300; $top = 230; $right = 2100; $bottom = 1250
  Draw-Axes $g $left $top $right $bottom "" "Diagnostic group" "Moran's I" 0 0.45
  $groups = @(
    @{ Scale="10m"; Label="10 m`nvariables"; Type="variable"; Color=$green },
    @{ Scale="20m"; Label="20 m`nvariables"; Type="variable"; Color=$green },
    @{ Scale="50m"; Label="50 m`nvariables"; Type="variable"; Color=$green },
    @{ Scale="10m"; Label="10 m`nresiduals"; Type="residual"; Color=$blue },
    @{ Scale="20m"; Label="20 m`nresiduals"; Type="residual"; Color=$blue },
    @{ Scale="50m"; Label="50 m`nresiduals"; Type="residual"; Color=$blue }
  )
  for ($i = 0; $i -lt $groups.Count; $i++) {
    $g1 = $groups[$i]
    $vals = @($moranKeep | Where-Object { $_.scale -eq $g1.Scale -and $_.diagnostic_type -eq $g1.Type } | ForEach-Object { [double]$_.moran_i })
    if ($vals.Count -eq 0) { continue }
    $avg = ($vals | Measure-Object -Average).Average
    $x = $left + 90 + $i * 280
    $y = Scale-Y $avg $top $bottom 0 0.45
    $g.FillRectangle((Brush ($g1.Color)), $x, $y, 150, $bottom - $y)
    Draw-CenteredText $g ("{0:0.000}" -f $avg) ($x - 30) ($y - 50) 210 44 24 $g1.Color ([System.Drawing.FontStyle]::Bold)
    Draw-CenteredText $g $g1.Label ($x - 45) ($bottom + 30) 240 80 22 $ink
  }
} | Out-Null

$transPriority = $trans | Where-Object { $_.x_metric -in @("phy_afaith","phy_rao","sp_shannon") -and $_.y_metric -in @("spec_spca_alpha","spec_spca_mean","spec_spca_rao") } | Sort-Object {[double]$_.delta_abs_r} -Descending | Select-Object -First 8
Save-Figure "08_transformation_sensitivity_priority_relationships.png" {
  param($g, $w, $h)
  Draw-Text $g "Transformation sensitivity for priority relationships" 100 70 42 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Text $g "Increase in absolute Pearson r relative to untransformed metrics" 100 130 27 $muted
  $left = 360; $top = 230; $right = 2120; $bottom = 1250
  Draw-Axes $g $left $top $right $bottom "" "Priority spectral-biodiversity relationship" "Delta absolute r" 0 0.16
  for ($i = 0; $i -lt $transPriority.Count; $i++) {
    $row = $transPriority[$i]
    $v = [double]$row.delta_abs_r
    $x = $left + 55 + $i * 205
    $y = Scale-Y $v $top $bottom 0 0.16
    $g.FillRectangle((Brush $purple), $x, $y, 108, $bottom - $y)
    Draw-CenteredText $g ("{0:0.000}" -f $v) ($x - 25) ($y - 50) 158 44 22 $purple ([System.Drawing.FontStyle]::Bold)
    $lab = ($row.scale.Replace("m"," m") + "`n" + (Label-Bio $row.x_metric).Replace("abundance-weighted ","A-wtd ").Replace("phylogenetic ","phylo. ") + "`n" + (Label-Spectral $row.y_metric).Replace("standardized PCA ",""))
    Draw-CenteredText $g $lab ($x - 62) ($bottom + 22) 230 145 18 $ink
  }
} | Out-Null

$tablePath = Join-Path $outDir "mdpi_remote_sensing_tables.md"
$mainRel = $all | Where-Object { $_.y_metric -in $specMetrics -and $_.x_metric -in @("phy_afaith","phy_rao","sp_shannon") } | Sort-Object scale, @{Expression={[math]::Abs([double]$_.pearson_r)}; Descending=$true}
$coverage = foreach ($scale in $scales) {
  $rows = @($all | Where-Object { $_.scale -eq $scale -and $_.y_metric -in $specMetrics -and $_.x_metric -in @("phy_afaith","phy_rao","sp_shannon") })
  [pscustomobject]@{
    Scale = $scale.Replace("m", " m")
    "Spectral-biodiversity complete cases (n)" = (($rows | ForEach-Object {[int]$_.n}) | Sort-Object -Unique) -join ", "
    "Primary spectral metrics" = "standardized PCA alpha hull; standardized PCA mean distance; standardized PCA Rao's Q"
    "Primary biodiversity contrast" = "abundance-weighted Faith's PD; phylogenetic Rao's Q; Shannon diversity"
  }
}
$topByScale = foreach ($scale in $scales) {
  $all | Where-Object { $_.scale -eq $scale -and $_.y_metric -in $specMetrics } |
    Sort-Object @{Expression={[math]::Abs([double]$_.pearson_r)}; Descending=$true} |
    Select-Object -First 6 |
    ForEach-Object {
      [pscustomobject]@{
        Scale = $scale.Replace("m"," m")
        "Biodiversity metric" = Label-Bio $_.x_metric
        "Spectral metric" = Label-Spectral $_.y_metric
        n = $_.n
        r = Fmt-R $_.pearson_r
        R2 = Fmt-R2 $_.r_squared
        p = Fmt-P $_.f_p_value
      }
    }
}
$bioConcord = $bioPairs | Where-Object {
  (($_.x_metric -eq "sp_shannon" -and $_.y_metric -eq "sp_simpson") -or
   ($_.x_metric -eq "sp_rich" -and $_.y_metric -eq "sp_shannon") -or
   ($_.x_metric -eq "phy_faith" -and $_.y_metric -eq "sp_rich") -or
   ($_.x_metric -eq "phy_rao" -and $_.y_metric -eq "phy_afaith") -or
   ($_.x_metric -eq "sp_shannon" -and $_.y_metric -eq "phy_afaith") -or
   ($_.x_metric -eq "phy_afaith" -and $_.y_metric -eq "sp_shannon"))
} | ForEach-Object {
  [pscustomobject]@{
    Scale = $_.scale.Replace("m"," m")
    "Metric pair" = "$(Label-Bio $_.x_metric) vs. $(Label-Bio $_.y_metric)"
    n = $_.n
    r = Fmt-R $_.pearson_r
    R2 = Fmt-R2 $_.r_squared
    p = Fmt-P $_.f_p_value
  }
}

function To-MdTable($rows) {
  if (-not $rows -or $rows.Count -eq 0) { return "_No rows available._`n" }
  $cols = @($rows[0].PSObject.Properties.Name)
  $out = @()
  $out += "| " + ($cols -join " | ") + " |"
  $out += "| " + (($cols | ForEach-Object { "---" }) -join " | ") + " |"
  foreach ($row in $rows) {
    $out += "| " + (($cols | ForEach-Object { [string]$row.$_ }) -join " | ") + " |"
  }
  return ($out -join "`r`n") + "`r`n"
}

$md = @()
$md += "# MDPI Remote Sensing Figure and Table Package"
$md += ""
$md += "Tables use manuscript-facing terminology from `SVH_Paper.docx`. Internal variable names are included only where they help connect the table to analysis outputs."
$md += ""
$md += "## Table 1. Metric definitions and calculation roles"
$md += ""
$md += To-MdTable @(
  [pscustomobject]@{ "Metric family"="Spectral"; "Manuscript label"="standardized PCA alpha hull"; "Variable"="spec_spca_alpha"; "Definition / role"="Area occupied by quadrat pixels in standardized PCA spectral space; emphasizes the spectral envelope." },
  [pscustomobject]@{ "Metric family"="Spectral"; "Manuscript label"="standardized PCA mean distance"; "Variable"="spec_spca_mean"; "Definition / role"="Mean Euclidean distance of quadrat pixels from the quadrat spectral centroid in standardized PCA space; emphasizes typical dispersion." },
  [pscustomobject]@{ "Metric family"="Spectral"; "Manuscript label"="standardized PCA Rao's Q"; "Variable"="spec_spca_rao"; "Definition / role"="Rao-style spectral dissimilarity in standardized PCA space; emphasizes pairwise squared spectral distances." },
  [pscustomobject]@{ "Metric family"="Biodiversity"; "Manuscript label"="abundance-weighted Faith's PD"; "Variable"="phy_afaith"; "Definition / role"="Phylogenetic diversity weighted by canopy crown-overlap abundance; primary phylogenetic response emphasized in the abstract." },
  [pscustomobject]@{ "Metric family"="Biodiversity"; "Manuscript label"="phylogenetic Rao's Q"; "Variable"="phy_rao"; "Definition / role"="Abundance-weighted phylogenetic dissimilarity using pairwise cophenetic distances among species." },
  [pscustomobject]@{ "Metric family"="Biodiversity"; "Manuscript label"="Shannon diversity"; "Variable"="sp_shannon"; "Definition / role"="Conventional species-diversity contrast based on proportional crown-overlap abundance." }
)
$md += ""
$md += "## Table 2. Quadrat coverage and complete-case counts"
$md += ""
$md += To-MdTable @($coverage)
$md += ""
$md += "## Table 3. Top standardized PCA spectral-biodiversity correlations by scale"
$md += ""
$md += To-MdTable @($topByScale)
$md += ""
$md += "## Table 4. Biodiversity metric concordance by scale"
$md += ""
$md += To-MdTable @($bioConcord)
$md += ""
$md += "## Table 5. Final main-versus-supplement figure inventory"
$md += ""
$md += To-MdTable @(
  [pscustomobject]@{ "Figure"="1"; "File"="01_study_site_quadrat_scales_workflow.png"; "Recommended placement"="Main"; "Purpose"="Introduces the site, shared 10, 20, and 50 m quadrat grains, and the analysis workflow." },
  [pscustomobject]@{ "Figure"="2"; "File"="02_standardized_pca_spectral_metric_relationships_by_scale.png"; "Recommended placement"="Main"; "Purpose"="Shows how standardized PCA spectral heterogeneity metrics relate to one another across scale." },
  [pscustomobject]@{ "Figure"="3"; "File"="03_biodiversity_metric_relationships_by_scale.png"; "Recommended placement"="Main"; "Purpose"="Shows concordance and divergence among species and phylogenetic biodiversity metrics." },
  [pscustomobject]@{ "Figure"="4"; "File"="04_standardized_pca_spectral_biodiversity_correlations_by_scale.png"; "Recommended placement"="Main"; "Purpose"="Summarizes all primary standardized PCA spectral-biodiversity correlations by scale." },
  [pscustomobject]@{ "Figure"="5"; "File"="05_abundance_weighted_faith_pd_scale_response.png"; "Recommended placement"="Main"; "Purpose"="Highlights the scale-dependent relationship with abundance-weighted Faith's PD." },
  [pscustomobject]@{ "Figure"="S1"; "File"="06_raw_vs_vector_normalized_pca_brightness_relationship.png"; "Recommended placement"="Supplement"; "Purpose"="Documents brightness dominance in raw PCA and reduction after vector normalization." },
  [pscustomobject]@{ "Figure"="S2"; "File"="07_spatial_moran_i_diagnostics.png"; "Recommended placement"="Supplement"; "Purpose"="Summarizes Moran's I diagnostics for priority variables and residuals." },
  [pscustomobject]@{ "Figure"="S3"; "File"="08_transformation_sensitivity_priority_relationships.png"; "Recommended placement"="Supplement"; "Purpose"="Shows that transformations mainly affected weaker priority relationships." }
)

[System.IO.File]::WriteAllText($tablePath, ($md -join "`r`n"), [System.Text.Encoding]::UTF8)
Write-Host "Created figures and tables in $outDir"
