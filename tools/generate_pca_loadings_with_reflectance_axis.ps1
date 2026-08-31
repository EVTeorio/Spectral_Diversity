Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$tableDir = Join-Path $root "reports/tables/pca_loading_spectral_regions"
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(22, 26, 32)
$muted = [System.Drawing.Color]::FromArgb(84, 91, 101)
$grid = [System.Drawing.Color]::FromArgb(218, 223, 230)
$pc1Color = [System.Drawing.Color]::FromArgb(24, 89, 153)
$pc2Color = [System.Drawing.Color]::FromArgb(185, 79, 43)
$reflectColor = [System.Drawing.Color]::FromArgb(112, 22, 22)

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

function Scale-X($x, $xmin, $xmax, $left, $right) {
  return $left + (([double]$x - [double]$xmin) / ([double]$xmax - [double]$xmin)) * ($right - $left)
}

function Scale-Y($y, $ymin, $ymax, $top, $bottom) {
  return $bottom - (([double]$y - [double]$ymin) / ([double]$ymax - [double]$ymin)) * ($bottom - $top)
}

function Get-Wavelength-Color($nm) {
  $w = [double]$nm
  if ($w -lt 450) { return [System.Drawing.Color]::FromArgb(92, 64, 172) }
  if ($w -lt 495) { return [System.Drawing.Color]::FromArgb(34, 96, 190) }
  if ($w -lt 570) { return [System.Drawing.Color]::FromArgb(36, 145, 77) }
  if ($w -lt 590) { return [System.Drawing.Color]::FromArgb(218, 181, 47) }
  if ($w -lt 620) { return [System.Drawing.Color]::FromArgb(226, 134, 39) }
  if ($w -lt 750) { return [System.Drawing.Color]::FromArgb(185, 41, 38) }
  return [System.Drawing.Color]::FromArgb(112, 22, 22)
}

function Add-Reflectance($rows, $reflectanceByWavelength) {
  foreach ($row in $rows) {
    $nmKey = [string]([int][double]$row.wavelength_nm)
    $row | Add-Member -NotePropertyName foliage_reflectance -NotePropertyValue ([double]$reflectanceByWavelength[$nmKey]) -Force
  }
  return $rows
}

function Draw-Series($g, $rows, $field, $xmin, $xmax, $ymin, $ymax, $left, $top, $right, $bottom, $color, $width = 4) {
  $pen = New-Pen $color $width
  for ($i = 1; $i -lt $rows.Count; $i++) {
    $x1 = Scale-X ([double]$rows[$i - 1].wavelength_nm) $xmin $xmax $left $right
    $y1 = Scale-Y ([double]$rows[$i - 1].$field) $ymin $ymax $top $bottom
    $x2 = Scale-X ([double]$rows[$i].wavelength_nm) $xmin $xmax $left $right
    $y2 = Scale-Y ([double]$rows[$i].$field) $ymin $ymax $top $bottom
    $g.DrawLine($pen, [single]$x1, [single]$y1, [single]$x2, [single]$y2)
  }
  $pen.Dispose()
}

function Draw-Reflectance-Series($g, $rows, $xmin, $xmax, $rmin, $rmax, $left, $top, $right, $bottom) {
  for ($i = 1; $i -lt $rows.Count; $i++) {
    $nm = ([double]$rows[$i - 1].wavelength_nm + [double]$rows[$i].wavelength_nm) / 2.0
    $pen = New-Pen (Get-Wavelength-Color $nm) 8
    $x1 = Scale-X ([double]$rows[$i - 1].wavelength_nm) $xmin $xmax $left $right
    $x2 = Scale-X ([double]$rows[$i].wavelength_nm) $xmin $xmax $left $right
    $y1 = Scale-Y ([double]$rows[$i - 1].foliage_reflectance) $rmin $rmax $top $bottom
    $y2 = Scale-Y ([double]$rows[$i].foliage_reflectance) $rmin $rmax $top $bottom
    $g.DrawLine($pen, [single]$x1, [single]$y1, [single]$x2, [single]$y2)
    $pen.Dispose()
  }
}

function Draw-Panel($g, $rows, $title, $left, $top, $right, $bottom, $loadingMin, $loadingMax, $reflectanceMin, $reflectanceMax) {
  $xmin = 398.0
  $xmax = 998.0

  $gridPen = New-Pen $grid 1
  $axisPen = New-Pen $ink 2
  $rightAxisPen = New-Pen $reflectColor 2
  $zeroPen = New-Pen ([System.Drawing.Color]::FromArgb(145, 88, 96, 105)) 2
  $zeroPen.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash

  foreach ($xTick in @(400, 500, 600, 700, 800, 900, 1000)) {
    $x = Scale-X $xTick $xmin $xmax $left $right
    $g.DrawLine($gridPen, [single]$x, [single]$top, [single]$x, [single]$bottom)
    Draw-Centered $g ("{0}" -f $xTick) ($x - 35) ($bottom + 10) 70 26 16 $muted
  }
  for ($yTick = -0.20; $yTick -le 0.20; $yTick += 0.10) {
    $y = Scale-Y $yTick $loadingMin $loadingMax $top $bottom
    $g.DrawLine($gridPen, [single]$left, [single]$y, [single]$right, [single]$y)
    Draw-Centered $g ("{0:0.00}" -f $yTick) ($left - 86) ($y - 13) 72 26 16 $muted
  }
  foreach ($rTick in @(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70)) {
    $y = Scale-Y $rTick $reflectanceMin $reflectanceMax $top $bottom
    Draw-Centered $g ("{0:0.00}" -f $rTick) ($right + 12) ($y - 13) 72 26 16 $reflectColor
  }

  $zeroY = Scale-Y 0 $loadingMin $loadingMax $top $bottom
  $g.DrawLine($zeroPen, [single]$left, [single]$zeroY, [single]$right, [single]$zeroY)

  Draw-Reflectance-Series $g $rows $xmin $xmax $reflectanceMin $reflectanceMax $left $top $right $bottom
  Draw-Series $g $rows "pc1_loading" $xmin $xmax $loadingMin $loadingMax $left $top $right $bottom $pc1Color 4.2
  Draw-Series $g $rows "pc2_loading" $xmin $xmax $loadingMin $loadingMax $left $top $right $bottom $pc2Color 4.2

  $g.DrawLine($axisPen, [single]$left, [single]$bottom, [single]$right, [single]$bottom)
  $g.DrawLine($axisPen, [single]$left, [single]$top, [single]$left, [single]$bottom)
  $g.DrawLine($rightAxisPen, [single]$right, [single]$top, [single]$right, [single]$bottom)

  Draw-Text $g $title ($left + 8) ($top - 48) 28 $ink ([System.Drawing.FontStyle]::Bold)
  $gridPen.Dispose()
  $axisPen.Dispose()
  $rightAxisPen.Dispose()
  $zeroPen.Dispose()
}

$loadings = @(Import-Csv -Path (Join-Path $tableDir "pca_pc1_pc2_loadings_by_wavelength.csv"))
$reflectanceRows = @(Import-Csv -Path (Join-Path $tableDir "sampled_plot_foliage_reflectance_50m.csv"))
$reflectanceByWavelength = @{}
foreach ($row in $reflectanceRows) {
  $reflectanceByWavelength[[string]([int][double]$row.wavelength_nm)] = [double]$row.mean_reflectance
}

$rawRows = Add-Reflectance @($loadings | Where-Object { $_.pca_basis -eq "regular_PCA" } | Sort-Object {[double]$_.wavelength_nm}) $reflectanceByWavelength
$stdRows = Add-Reflectance @($loadings | Where-Object { $_.pca_basis -eq "standardized_PCA" } | Sort-Object {[double]$_.wavelength_nm}) $reflectanceByWavelength

$allLoadings = @()
foreach ($row in @($rawRows + $stdRows)) {
  $allLoadings += [double]$row.pc1_loading
  $allLoadings += [double]$row.pc2_loading
}
$loadingAbs = (($allLoadings | ForEach-Object { [Math]::Abs($_) }) | Measure-Object -Maximum).Maximum
$loadingLimit = [Math]::Ceiling(($loadingAbs + 0.01) * 10) / 10
$loadingMin = -1 * $loadingLimit
$loadingMax = $loadingLimit

$allReflectance = @($reflectanceRows | ForEach-Object { [double]$_.mean_reflectance })
$reflectanceMin = 0.0
$reflectanceMax = 0.70

$width = 2450
$height = 1850
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "PCA loadings across wavelength" 145 80 42 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "PC1 and PC2 loadings with plot-level foliage reflectance shown on the right axis" 145 140 24 $muted

$left = 250
$right = 2195
$top1 = 330
$bottom1 = 840
$top2 = 1120
$bottom2 = 1630

Draw-Panel $g $rawRows "Raw PCA" $left $top1 $right $bottom1 $loadingMin $loadingMax $reflectanceMin $reflectanceMax
Draw-Panel $g $stdRows "Vector-normalized PCA" $left $top2 $right $bottom2 $loadingMin $loadingMax $reflectanceMin $reflectanceMax

Draw-Centered $g "Wavelength (nm)" 775 1740 900 44 27 $ink ([System.Drawing.FontStyle]::Bold)

$font = New-Font 27 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $ink
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("PCA loading", $font, $brush, [System.Drawing.RectangleF]::new(55, 680, 54, 600), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

$font = New-Font 24 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $reflectColor
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("Foliage reflectance", $font, $brush, [System.Drawing.RectangleF]::new(2330, 630, 54, 650), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

$legendX = 1410
$legendY = 95
$pc1Pen = New-Pen $pc1Color 5
$pc2Pen = New-Pen $pc2Color 5
$reflectPen = New-Pen $reflectColor 8

$g.DrawLine($pc1Pen, [single]$legendX, [single]$legendY, [single]($legendX + 90), [single]$legendY)
Draw-Text $g "PC1 loading" ($legendX + 112) ($legendY - 20) 20 $ink
$g.DrawLine($pc2Pen, [single]$legendX, [single]($legendY + 45), [single]($legendX + 90), [single]($legendY + 45))
Draw-Text $g "PC2 loading" ($legendX + 112) ($legendY + 25) 20 $ink
$g.DrawLine($reflectPen, [single]$legendX, [single]($legendY + 90), [single]($legendX + 90), [single]($legendY + 90))
Draw-Text $g "Plot foliage reflectance" ($legendX + 112) ($legendY + 70) 20 $ink

$pc1Pen.Dispose()
$pc2Pen.Dispose()
$reflectPen.Dispose()

$outPath = Join-Path $outDir "17_pca_loadings_with_foliage_reflectance_axis.png"
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Output $outPath
