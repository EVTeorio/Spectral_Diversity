Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$tile = "sub50_16"
$hsiPath = Join-Path $root "Quad_Spectra/50m_smooth_5nm/$tile"
$outDir = Join-Path $root "Documents/Tables and Figures"
$outPath = Join-Path $outDir "18_rgb_shadow_mask_before_after_50m_tile.png"

New-Item -ItemType Directory -Force -Path $outDir | Out-Null
Add-Type -AssemblyName System.Drawing

$samples = 680
$lines = 680
$n = $samples * $lines
$threshold = 0.0305476

function Band-Index([int]$wavelength) {
  return [int](($wavelength - 398) / 5) + 1
}

function Read-EnviBand([string]$path, [int]$bandIndex1) {
  $offset = [int64]($bandIndex1 - 1) * $script:n * 4
  $vals = [single[]]::new($script:n)
  $fs = [System.IO.File]::Open($path, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)
  try {
    [void]$fs.Seek($offset, [System.IO.SeekOrigin]::Begin)
    $br = [System.IO.BinaryReader]::new($fs)
    for ($i = 0; $i -lt $script:n; $i++) {
      $vals[$i] = $br.ReadSingle()
    }
  }
  finally {
    $fs.Close()
  }
  return $vals
}

function Percentile([single[]]$values, [double]$p) {
  $valid = [System.Collections.Generic.List[double]]::new()
  foreach ($v in $values) {
    if (-not [double]::IsNaN($v)) {
      $valid.Add([double]$v)
    }
  }
  $arr = $valid.ToArray()
  [Array]::Sort($arr)
  if ($arr.Length -eq 0) { return 0.0 }
  $idx = [int][Math]::Round(($arr.Length - 1) * $p)
  return $arr[[Math]::Max(0, [Math]::Min($arr.Length - 1, $idx))]
}

function Scale-Byte([double]$v, [double]$lo, [double]$hi) {
  if ([double]::IsNaN($v)) { return 0 }
  if ($hi -le $lo) { return 0 }
  $x = 255.0 * (($v - $lo) / ($hi - $lo))
  if ($x -lt 0) { return 0 }
  if ($x -gt 255) { return 255 }
  return [byte][Math]::Round($x)
}

function New-RgbBitmap([single[]]$red, [single[]]$green, [single[]]$blue, [single[]]$maskBand) {
  $rLo = Percentile $red 0.02
  $rHi = Percentile $red 0.98
  $gLo = Percentile $green 0.02
  $gHi = Percentile $green 0.98
  $bLo = Percentile $blue 0.02
  $bHi = Percentile $blue 0.98

  $bmp = [System.Drawing.Bitmap]::new($script:samples, $script:lines, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  $rect = [System.Drawing.Rectangle]::new(0, 0, $script:samples, $script:lines)
  $data = $bmp.LockBits($rect, [System.Drawing.Imaging.ImageLockMode]::WriteOnly, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  try {
    $bytes = [byte[]]::new($data.Stride * $script:lines)
    for ($y = 0; $y -lt $script:lines; $y++) {
      for ($x = 0; $x -lt $script:samples; $x++) {
        $i = $y * $script:samples + $x
        $j = $y * $data.Stride + $x * 4
        $valid = (-not [double]::IsNaN($red[$i])) -and (-not [double]::IsNaN($green[$i])) -and (-not [double]::IsNaN($blue[$i])) -and (-not [double]::IsNaN($maskBand[$i]))
        $illuminated = $valid -and ([double]$maskBand[$i] -gt $script:threshold)
        if ($valid) {
          $bytes[$j + 0] = Scale-Byte $blue[$i] $bLo $bHi
          $bytes[$j + 1] = Scale-Byte $green[$i] $gLo $gHi
          $bytes[$j + 2] = Scale-Byte $red[$i] $rLo $rHi
          $bytes[$j + 3] = 255
        }
      }
    }
    [System.Runtime.InteropServices.Marshal]::Copy($bytes, 0, $data.Scan0, $bytes.Length)
  }
  finally {
    $bmp.UnlockBits($data)
  }
  return $bmp
}

function New-HatchOverlay([single[]]$maskBand) {
  $bmp = [System.Drawing.Bitmap]::new($script:samples, $script:lines, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  $rect = [System.Drawing.Rectangle]::new(0, 0, $script:samples, $script:lines)
  $data = $bmp.LockBits($rect, [System.Drawing.Imaging.ImageLockMode]::WriteOnly, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
  try {
    $bytes = [byte[]]::new($data.Stride * $script:lines)
    for ($y = 0; $y -lt $script:lines; $y++) {
      for ($x = 0; $x -lt $script:samples; $x++) {
        $i = $y * $script:samples + $x
        $j = $y * $data.Stride + $x * 4
        $valid = -not [double]::IsNaN($maskBand[$i])
        $shadow = $valid -and ([double]$maskBand[$i] -le $script:threshold)
        if ($shadow) {
          $bytes[$j + 0] = 205
          $bytes[$j + 1] = 70
          $bytes[$j + 2] = 170
          $bytes[$j + 3] = 132
          if ((($x + $y) % 16) -lt 4) {
            $bytes[$j + 0] = 28
            $bytes[$j + 1] = 32
            $bytes[$j + 2] = 238
            $bytes[$j + 3] = 255
          }
        }
      }
    }
    [System.Runtime.InteropServices.Marshal]::Copy($bytes, 0, $data.Scan0, $bytes.Length)
  }
  finally {
    $bmp.UnlockBits($data)
  }
  return $bmp
}

function Font([single]$size, $style = [System.Drawing.FontStyle]::Regular) {
  return [System.Drawing.Font]::new("Arial", $size, $style)
}

function Brush($color) {
  return [System.Drawing.SolidBrush]::new($color)
}

function PenC($color, [single]$width = 2) {
  return [System.Drawing.Pen]::new($color, $width)
}

function Draw-Text($g, [string]$text, [single]$x, [single]$y, [single]$size, $color, $style = [System.Drawing.FontStyle]::Regular) {
  $font = Font $size $style
  $brush = Brush $color
  $g.DrawString($text, $font, $brush, $x, $y)
  $font.Dispose()
  $brush.Dispose()
}

function Draw-Centered($g, [string]$text, [single]$x, [single]$y, [single]$w, [single]$h, [single]$size, $color, $style = [System.Drawing.FontStyle]::Regular) {
  $font = Font $size $style
  $brush = Brush $color
  $fmt = [System.Drawing.StringFormat]::new()
  $fmt.Alignment = [System.Drawing.StringAlignment]::Center
  $fmt.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $font, $brush, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $fmt)
  $fmt.Dispose()
  $font.Dispose()
  $brush.Dispose()
}

$red = Read-EnviBand $hsiPath (Band-Index 648)
$green = Read-EnviBand $hsiPath (Band-Index 548)
$blue = Read-EnviBand $hsiPath (Band-Index 453)
$maskBand = Read-EnviBand $hsiPath (Band-Index 563)

$valid = 0
$illuminated = 0
foreach ($v in $maskBand) {
  if (-not [double]::IsNaN($v)) {
    $valid++
    if ([double]$v -gt $threshold) { $illuminated++ }
  }
}
$shadow = $valid - $illuminated
$illumPct = 100.0 * $illuminated / [Math]::Max(1, $valid)
$shadowPct = 100.0 * $shadow / [Math]::Max(1, $valid)

$before = New-RgbBitmap $red $green $blue $maskBand
$after = New-RgbBitmap $red $green $blue $maskBand
$hatch = New-HatchOverlay $maskBand

$w = 2400
$h = 1320
$fig = [System.Drawing.Bitmap]::new($w, $h, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($fig)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

$ink = [System.Drawing.Color]::FromArgb(38, 45, 52)
$muted = [System.Drawing.Color]::FromArgb(96, 104, 112)
$blueInk = [System.Drawing.Color]::FromArgb(44, 92, 152)
$greenInk = [System.Drawing.Color]::FromArgb(38, 122, 80)
$shadowFill = [System.Drawing.Color]::FromArgb(132, 170, 70, 205)
$shadowInk = [System.Drawing.Color]::FromArgb(255, 238, 32, 28)
$border = [System.Drawing.Color]::FromArgb(170, 178, 186)

$g.CompositingMode = [System.Drawing.Drawing2D.CompositingMode]::SourceOver
Draw-Text $g "RGB shadow-masking example for a 50 m hyperspectral tile" 90 55 40 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Illuminated pixels retained where 563 nm reflectance > 0.0305476" 90 112 27 $muted

$panelTop = 235
$panelSize = 850
$leftX = 135
$rightX = 1405
$labelY = 178

Draw-Centered $g "A. RGB tile before masking" $leftX $labelY $panelSize 48 29 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Centered $g "B. RGB tile with shadow mask" $rightX $labelY $panelSize 48 29 $ink ([System.Drawing.FontStyle]::Bold)

$g.DrawImage($before, [System.Drawing.Rectangle]::new($leftX, $panelTop, $panelSize, $panelSize))
$g.DrawImage($after, [System.Drawing.Rectangle]::new($rightX, $panelTop, $panelSize, $panelSize))
$g.DrawImage($hatch, [System.Drawing.Rectangle]::new($rightX, $panelTop, $panelSize, $panelSize))
$g.DrawRectangle((PenC $border 3), $leftX, $panelTop, $panelSize, $panelSize)
$g.DrawRectangle((PenC $border 3), $rightX, $panelTop, $panelSize, $panelSize)

$arrowPen = PenC $muted 7
$arrowPen.EndCap = [System.Drawing.Drawing2D.LineCap]::ArrowAnchor
$g.DrawLine($arrowPen, $leftX + $panelSize + 90, $panelTop + 430, $rightX - 90, $panelTop + 430)
$arrowPen.Dispose()
Draw-Centered $g "shadow mask" ($leftX + $panelSize + 118) ($panelTop + 342) 250 60 25 $muted ([System.Drawing.FontStyle]::Bold)

$legendX = 185
$legendY = 1140
$sw = 42
$g.FillRectangle((Brush $greenInk), $legendX, $legendY, $sw, $sw)
Draw-Text $g ("Illuminated pixels retained: {0:N1}%" -f $illumPct) ($legendX + 62) ($legendY - 2) 25 $ink
$legendX2 = 1020
$g.FillRectangle((Brush $shadowFill), $legendX2, $legendY, $sw, $sw)
$checkerPen = PenC $shadowInk 5
$clipState = $g.Save()
$g.SetClip([System.Drawing.Rectangle]::new($legendX2 + 2, $legendY + 2, $sw - 4, $sw - 4))
for ($d = -$sw; $d -le $sw; $d += 13) {
  $g.DrawLine($checkerPen, $legendX2 + $d, $legendY + $sw, $legendX2 + $d + $sw, $legendY)
}
$g.Restore($clipState)
$g.DrawRectangle((PenC $shadowInk 3), $legendX2, $legendY, $sw, $sw)
$checkerPen.Dispose()
Draw-Text $g ("Shadowed pixels masked: {0:N1}%" -f $shadowPct) ($legendX2 + 62) ($legendY - 2) 25 $ink

$fig.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)

$g.Dispose()
$fig.Dispose()
$before.Dispose()
$after.Dispose()
$hatch.Dispose()

Write-Host "Created $outPath"
Write-Host ("Valid pixels: {0}; retained illuminated: {1} ({2:N1}%); masked shadow: {3} ({4:N1}%)" -f $valid, $illuminated, $illumPct, $shadow, $shadowPct)
