Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(22, 26, 32)
$muted = [System.Drawing.Color]::FromArgb(86, 94, 104)
$grid20 = [System.Drawing.Color]::FromArgb(245, 54, 72, 88)
$grid10 = [System.Drawing.Color]::FromArgb(235, 34, 114, 164)
$grid50 = [System.Drawing.Color]::FromArgb(255, 18, 22, 28)
$callout = [System.Drawing.Color]::FromArgb(220, 255, 255, 255)

function New-Font($size, $style = [System.Drawing.FontStyle]::Regular) {
  return [System.Drawing.Font]::new($fontFamily, [single]$size, $style)
}

function New-Brush($color) {
  return [System.Drawing.SolidBrush]::new($color)
}

function New-Pen($color, $width = 1) {
  return [System.Drawing.Pen]::new($color, [single]$width)
}

function Draw-Text($g, $text, $x, $y, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $g.DrawString($text, $font, $brush, [single]$x, [single]$y)
  $font.Dispose()
  $brush.Dispose()
}

function Draw-Centered($g, $text, $x, $y, $w, $h, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $fmt = [System.Drawing.StringFormat]::new()
  $fmt.Alignment = [System.Drawing.StringAlignment]::Center
  $fmt.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $font, $brush, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $fmt)
  $fmt.Dispose()
  $font.Dispose()
  $brush.Dispose()
}

function Draw-LegendSpecies($g, $name, $count, $x, $y, $size = 20) {
  $italicFont = New-Font $size ([System.Drawing.FontStyle]::Italic)
  $regularFont = New-Font $size ([System.Drawing.FontStyle]::Regular)
  $brush = New-Brush $ink
  $g.DrawString($name, $italicFont, $brush, [single]$x, [single]$y)
  $nameWidth = $g.MeasureString($name, $italicFont).Width
  $g.DrawString((" (" + $count + ")"), $regularFont, $brush, [single]($x + $nameWidth - 2), [single]$y)
  $italicFont.Dispose()
  $regularFont.Dispose()
  $brush.Dispose()
}

function Get-ScientificName($taxonomy, $code) {
  if ($taxonomy.ContainsKey($code)) {
    return $taxonomy[$code]
  }
  return $code
}

function Get-SpeciesColor($code) {
  $hash = 0
  foreach ($ch in $code.ToCharArray()) {
    $hash = (($hash * 31) + [int][char]$ch) % 360
  }
  $hue = [double]$hash
  $sat = 0.58
  $val = 0.82
  $c = $val * $sat
  $x = $c * (1 - [Math]::Abs((($hue / 60.0) % 2) - 1))
  $m = $val - $c
  $r1 = 0.0
  $g1 = 0.0
  $b1 = 0.0
  if ($hue -lt 60) { $r1 = $c; $g1 = $x; $b1 = 0 }
  elseif ($hue -lt 120) { $r1 = $x; $g1 = $c; $b1 = 0 }
  elseif ($hue -lt 180) { $r1 = 0; $g1 = $c; $b1 = $x }
  elseif ($hue -lt 240) { $r1 = 0; $g1 = $x; $b1 = $c }
  elseif ($hue -lt 300) { $r1 = $x; $g1 = 0; $b1 = $c }
  else { $r1 = $c; $g1 = 0; $b1 = $x }
  return [System.Drawing.Color]::FromArgb(
    175,
    [int](255 * ($r1 + $m)),
    [int](255 * ($g1 + $m)),
    [int](255 * ($b1 + $m))
  )
}

function To-X($x, $x0, $scale, $left) {
  return $left + (([double]$x - [double]$x0) * $scale)
}

function To-Y($y, $y0, $scale, $bottom) {
  return $bottom - (([double]$y - [double]$y0) * $scale)
}

$x0 = 0.0
$y0 = 0.0
$quadratSize = 50.0
$treePath = Join-Path $root "PR_tree_DL.csv"
$taxaPath = Join-Path $root "51sp_taxanomy.csv"

$taxonomy = @{}
Import-Csv -Path $taxaPath | ForEach-Object {
  $speciesName = $_.species
  if ($speciesName -eq "sp") {
    $taxonomy[$_.sp_code] = $_.genus + " sp."
  } else {
    $taxonomy[$_.sp_code] = $_.genus + " " + $speciesName
  }
}

$trees = @(Import-Csv -Path $treePath | Where-Object {
  ($_.'DBH.2024' -as [double]) -ne $null -and
  [double]$_.'DBH.2024' -ge 200 -and
  $_.'crown.position' -in @("3", "4", "5") -and
  $_.cluster_status -in @("A", "R") -and
  $_.cw_m_2025 -ne "" -and $_.cw_m_2025 -ne "NA" -and
  ($_.'new.gx' -as [double]) -ne $null -and
  ($_.'new.gy' -as [double]) -ne $null -and
  [double]$_.cw_m_2025 -gt 0
} | ForEach-Object {
  $cx = [double]$_.'new.gx'
  $cy = [double]$_.'new.gy'
  $radius = [double]$_.cw_m_2025 / 2.0
  if (($cx + $radius) -ge $x0 -and ($cx - $radius) -le ($x0 + $quadratSize) -and
      ($cy + $radius) -ge $y0 -and ($cy - $radius) -le ($y0 + $quadratSize)) {
    [pscustomobject]@{
      Species = $_.sp
      X = $cx
      Y = $cy
      Radius = $radius
      DBH = [double]$_.'DBH.2024'
    }
  }
})

$width = 2400
$height = 1900
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "Tree crown intersections across nested quadrat grains" 145 80 42 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Southwest 50m quadrat" 145 140 24 $muted

$plotLeft = 335
$plotTop = 275
$plotSize = 1325
$plotBottom = $plotTop + $plotSize
$scale = $plotSize / $quadratSize

$clipPath = [System.Drawing.Drawing2D.GraphicsPath]::new()
$clipPath.AddRectangle([System.Drawing.RectangleF]::new($plotLeft, $plotTop, $plotSize, $plotSize))
$oldClip = $g.Clip
$g.SetClip($clipPath)

foreach ($tree in ($trees | Sort-Object Radius -Descending)) {
  $color = Get-SpeciesColor $tree.Species
  $brush = New-Brush $color
  $penColor = [System.Drawing.Color]::FromArgb(165, [int]($color.R * 0.62), [int]($color.G * 0.62), [int]($color.B * 0.62))
  $pen = New-Pen $penColor 1.2
  $cx = To-X $tree.X $x0 $scale $plotLeft
  $cy = To-Y $tree.Y $y0 $scale $plotBottom
  $r = $tree.Radius * $scale
  $g.FillEllipse($brush, [single]($cx - $r), [single]($cy - $r), [single](2 * $r), [single](2 * $r))
  $g.DrawEllipse($pen, [single]($cx - $r), [single]($cy - $r), [single](2 * $r), [single](2 * $r))
  $brush.Dispose()
  $pen.Dispose()
}

$g.Clip = $oldClip
$clipPath.Dispose()

$pen50 = New-Pen $grid50 6
$pen20 = New-Pen $grid20 4
$pen10 = New-Pen $grid10 3
$pen20.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Solid
$pen10.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash

$g.DrawRectangle($pen50, [single]$plotLeft, [single]$plotTop, [single]$plotSize, [single]$plotSize)

foreach ($offset in @(20, 40)) {
  $x = To-X ($x0 + $offset) $x0 $scale $plotLeft
  $g.DrawLine($pen20, [single]$x, [single](To-Y ($y0 + 0) $y0 $scale $plotBottom), [single]$x, [single](To-Y ($y0 + 40) $y0 $scale $plotBottom))
  $y = To-Y ($y0 + $offset) $y0 $scale $plotBottom
  $g.DrawLine($pen20, [single](To-X ($x0 + 0) $x0 $scale $plotLeft), [single]$y, [single](To-X ($x0 + 40) $x0 $scale $plotLeft), [single]$y)
}

$x10 = To-X ($x0 + 10) $x0 $scale $plotLeft
$y10 = To-Y ($y0 + 10) $y0 $scale $plotBottom
$g.DrawLine($pen10, [single]$x10, [single](To-Y ($y0 + 0) $y0 $scale $plotBottom), [single]$x10, [single](To-Y ($y0 + 20) $y0 $scale $plotBottom))
$g.DrawLine($pen10, [single](To-X ($x0 + 0) $x0 $scale $plotLeft), [single]$y10, [single](To-X ($x0 + 20) $x0 $scale $plotLeft), [single]$y10)

$calloutBrush = New-Brush $callout
$g.FillRectangle($calloutBrush, [single]($plotLeft + 22), [single]($plotTop + 22), 368, 100)
Draw-Text $g "50 m quadrat" ($plotLeft + 42) ($plotTop + 35) 28 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g ("Intersecting crowns: " + $trees.Count) ($plotLeft + 42) ($plotTop + 76) 21 $muted

$g.FillRectangle($calloutBrush, [single](To-X 1 $x0 $scale $plotLeft), [single](To-Y 39 $y0 $scale $plotBottom), 345, 76)
Draw-Text $g "Four 20 m quadrats" (To-X 2 $x0 $scale $plotLeft) ((To-Y 39 $y0 $scale $plotBottom) + 10) 22 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "aligned to plot edge" (To-X 2 $x0 $scale $plotLeft) ((To-Y 39 $y0 $scale $plotBottom) + 42) 18 $muted

$g.FillRectangle($calloutBrush, [single](To-X 21 $x0 $scale $plotLeft), [single](To-Y 18 $y0 $scale $plotBottom), 315, 76)
Draw-Text $g "One 20 m quadrat" (To-X 22 $x0 $scale $plotLeft) ((To-Y 18 $y0 $scale $plotBottom) + 10) 22 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "subdivided into 10 m" (To-X 22 $x0 $scale $plotLeft) ((To-Y 18 $y0 $scale $plotBottom) + 42) 18 $muted

$legendX = 1745
$legendY = 350
Draw-Text $g "Line key" $legendX ($legendY - 70) 28 $ink ([System.Drawing.FontStyle]::Bold)
$g.DrawLine($pen50, [single]$legendX, [single]$legendY, [single]($legendX + 125), [single]$legendY)
Draw-Text $g "50 m quadrat boundary" ($legendX + 155) ($legendY - 19) 22 $ink
$g.DrawLine($pen20, [single]$legendX, [single]($legendY + 70), [single]($legendX + 125), [single]($legendY + 70))
Draw-Text $g "20 m quadrat boundaries" ($legendX + 155) ($legendY + 51) 22 $ink
$g.DrawLine($pen10, [single]$legendX, [single]($legendY + 140), [single]($legendX + 125), [single]($legendY + 140))
Draw-Text $g "10 m subdivisions" ($legendX + 155) ($legendY + 121) 22 $ink

Draw-Text $g "Crown color = species" $legendX ($legendY + 255) 28 $ink ([System.Drawing.FontStyle]::Bold)
$speciesCounts = $trees | Group-Object Species | Sort-Object Count -Descending
$legendRows = @($speciesCounts | Select-Object -First 9)
$otherCount = 0
if ($speciesCounts.Count -gt 9) {
  $otherCount = ($speciesCounts | Select-Object -Skip 9 | Measure-Object Count -Sum).Sum
}

$rowY = $legendY + 310
foreach ($item in $legendRows) {
  $swatchColor = Get-SpeciesColor $item.Name
  $brush = New-Brush $swatchColor
  $g.FillRectangle($brush, [single]$legendX, [single]$rowY, 36, 24)
  $brush.Dispose()
  Draw-LegendSpecies $g (Get-ScientificName $taxonomy $item.Name) $item.Count ($legendX + 55) ($rowY - 7) 20
  $rowY += 42
}
if ($otherCount -gt 0) {
  Draw-Text $g ("Other species (" + $otherCount + ")") ($legendX + 55) ($rowY - 7) 20 $ink
}

Draw-Centered $g "Local plot x coordinate (m)" ($plotLeft + 215) ($plotBottom + 85) 895 45 25 $ink ([System.Drawing.FontStyle]::Bold)
$font = New-Font 25 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $ink
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("Local plot y coordinate (m)", $font, $brush, [System.Drawing.RectangleF]::new(110, $plotTop + 270, 54, 780), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

for ($i = 0; $i -le 5; $i++) {
  $tick = $i * 10
  $x = To-X $tick $x0 $scale $plotLeft
  $g.DrawLine($pen50, [single]$x, [single]$plotBottom, [single]$x, [single]($plotBottom + 10))
  Draw-Centered $g ("{0}" -f $tick) ($x - 35) ($plotBottom + 18) 70 26 17 $muted
  $y = To-Y $tick $y0 $scale $plotBottom
  $g.DrawLine($pen50, [single]($plotLeft - 10), [single]$y, [single]$plotLeft, [single]$y)
  Draw-Centered $g ("{0}" -f $tick) ($plotLeft - 80) ($y - 13) 58 26 17 $muted
}

$pen50.Dispose()
$pen20.Dispose()
$pen10.Dispose()
$calloutBrush.Dispose()

$outPath = Join-Path $outDir "12_methods_crown_quadrat_intersection_nested_scales.png"
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Output $outPath
